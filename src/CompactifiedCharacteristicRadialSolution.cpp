#include "CompactifiedCharacteristicRadialSolution.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {
using State = std::array<complex_type,2>;
float_type minimum(const float_type& x,const float_type& y) { return x<y?x:y; }
float_type maximum(const float_type& x,const float_type& y) { return x>y?x:y; }
bool finite_complex(const complex_type& x) {
    return mp::isfinite(x.real())&&mp::isfinite(x.imag());
}
std::vector<complex_type> product(const std::vector<complex_type>& a,
                                const std::vector<complex_type>& b) {
    std::vector<complex_type> c(a.size()+b.size()-1,complex_type(0));
    for(std::size_t k=0;k<a.size();++k)
        for(std::size_t j=0;j<b.size();++j) c[k+j]+=a[k]*b[j];
    return c;
}
State polynomial(const std::vector<complex_type>& c,const float_type& x) {
    State v{{c.back(),complex_type(0)}};
    for(std::size_t k=c.size()-1;k>0;--k) {
        v[1]=v[1]*x+v[0]; v[0]=v[0]*x+c[k-1];
    }
    return v;
}
complex_type polynomial_value(const std::vector<complex_type>& c,const float_type& x) {
    complex_type value=c.back();
    for(std::size_t n=c.size()-1;n>0;--n)value=value*x+c[n-1];
    return value;
}
std::array<complex_type,3> polynomial_derivatives(
    const std::vector<complex_type>& c,const float_type& x) {
    std::array<complex_type,3> value{{c.back(),0,0}};
    for(std::size_t n=c.size()-1;n>0;--n) {
        value[2]=value[2]*x+2*value[1];
        value[1]=value[1]*x+value[0];value[0]=value[0]*x+c[n-1];
    }
    return value;
}
float_type state_size(const State& x,const float_type& scale) {
    return maximum(float_type(mp::abs(x[0])),float_type(mp::abs(x[1]))*scale);
}
float_type term_size(const complex_type& c,unsigned n,const float_type& power) {
    return mp::abs(c)*power*float_type(std::max(1u,n));
}
}

CompactifiedCharacteristicRadialSolution::CompactifiedCharacteristicRadialSolution(
    float_type a,float_type lambda,float_type omega,int l,int m,float_type lo,float_type hi)
    : CompactifiedCharacteristicRadialSolution(a,lambda,omega,l,m,lo,hi,Options{}) {}

CompactifiedCharacteristicRadialSolution::CompactifiedCharacteristicRadialSolution(
    float_type a,float_type lambda,float_type omega,int l,int m,float_type lo,float_type hi,
    const Options& options)
    : a_(a),lambda_(lambda),omega_(omega),r_min_(lo),r_max_(hi),l_(l),m_(m),
      horizon_reference_(0),options_(options) {
    validate(); initialize_endpoints();
}

void CompactifiedCharacteristicRadialSolution::validate() {
    if(!mp::isfinite(a_)||!mp::isfinite(lambda_)||!mp::isfinite(omega_)||
       !mp::isfinite(r_min_)||!mp::isfinite(r_max_)||mp::abs(a_)>=1)
        throw std::invalid_argument("compactified scalar backend requires finite inputs and |a|<1");
    if(l_<0||m_>l_||m_<-l_)
        throw std::invalid_argument("compactified scalar backend requires l>=0 and |m|<=l");
    float_type root=mp::sqrt(1-a_*a_);
    rp_=1+root; rm_=1-root; b_=rm_/rp_; j_=a_*m_/rp_; w_=omega_*rp_;
    if(r_min_<=rp_||r_max_<r_min_)
        throw std::invalid_argument("compactified source interval must be outside the horizon");
    if(!mp::isfinite(options_.endpoint_tolerance)||
       !mp::isfinite(options_.integration_absolute_tolerance)||
       !mp::isfinite(options_.integration_relative_tolerance)||
       !mp::isfinite(options_.interpolation_tolerance)||
       !mp::isfinite(options_.outgoing_tolerance)||
       options_.endpoint_tolerance<=0||options_.integration_absolute_tolerance<=0||
       options_.integration_relative_tolerance<=0||options_.interpolation_tolerance<=0||
       options_.outgoing_tolerance<0||options_.horizon_order<8||
       options_.infinity_max_order<12||options_.max_steps==0)
        throw std::invalid_argument("invalid compactified startup or integration options");
    if(omega_==0 && lambda_!=float_type(l_)*(float_type(l_)+1))
        throw std::invalid_argument("static scalar lambda must equal l(l+1) exactly");
    sigma_min_=rp_/r_max_; sigma_max_=rp_/r_min_;
    horizon_reference_=HorizonStaticReference(l_,b_,j_);
}

float_type CompactifiedCharacteristicRadialSolution::tortoise(const float_type& a,
                                                               const float_type& r) {
    if(!mp::isfinite(a)||!mp::isfinite(r)||mp::abs(a)>=1)
        throw std::invalid_argument("tortoise requires finite subextremal exterior inputs");
    float_type root=mp::sqrt(1-a*a),rp=1+root,rm=1-root;
    if(root<=0||r<=rp) throw std::invalid_argument("tortoise requires subextremal exterior");
    // STORM normalizes the logarithms by 2.
    return r+(rp*mp::log((r-rp)/2)-rm*mp::log((r-rm)/2))/root;
}
float_type CompactifiedCharacteristicRadialSolution::azimuthal_phase(const float_type& a,
                                                                      const float_type& r) {
    if(!mp::isfinite(a)||!mp::isfinite(r)||mp::abs(a)>=1)
        throw std::invalid_argument("chi requires finite subextremal exterior inputs");
    float_type root=mp::sqrt(1-a*a),rp=1+root;
    if(root<=0||r<=rp) throw std::invalid_argument("chi requires subextremal exterior");
    float_type log_ratio=r-rp<root?float_type(mp::log(r-rp)-mp::log(r-(1-root))):
        float_type(mp::log1p(-2*root/(r-(1-root))));
    return a/(2*root)*log_ratio;
}

void CompactifiedCharacteristicRadialSolution::operator_coefficients(
    const float_type& a,const float_type& lambda,const float_type& omega,int m,int H,
    const float_type& sigma,complex_type& A,complex_type& B,complex_type& C) {
    if(H!=-1&&H!=1) throw std::invalid_argument("characteristic H must be -1 or +1");
    if(!mp::isfinite(a)||!mp::isfinite(lambda)||!mp::isfinite(omega)||
       !mp::isfinite(sigma)||mp::abs(a)>=1||sigma<0||sigma>1)
        throw std::invalid_argument("invalid compactified coefficient arguments");
    float_type root=mp::sqrt(1-a*a),rp=1+root,b=(1-root)/rp,j=a*m/rp,w=omega*rp;
    float_type s2=sigma*sigma,D=(1-sigma)*(1-b*sigma);
    A=s2*D;
    B=2*sigma-3*(1+b)*s2+4*b*s2*sigma-
        2*constants::I*H*(w*(1+b*s2)-j*s2);
    C=-lambda-(1+b)*sigma+2*b*s2-2*constants::I*H*(b*w-j)*sigma;
}
void CompactifiedCharacteristicRadialSolution::coefficients(int H,const float_type& sigma,
    complex_type& A,complex_type& B,complex_type& C) const {
    if(H!=-1&&H!=1)throw std::invalid_argument("characteristic H must be -1 or +1");
    if(!mp::isfinite(sigma)||sigma<0||sigma>1)
        throw std::invalid_argument("invalid compactified coefficient sigma");
    float_type s2=sigma*sigma,D=(1-sigma)*(1-b_*sigma);
    A=s2*D;
    B=2*sigma-3*(1+b_)*s2+4*b_*s2*sigma-
        2*constants::I*H*(w_*(1+b_*s2)-j_*s2);
    C=-lambda_-(1+b_)*sigma+2*b_*s2-2*constants::I*H*(b_*w_-j_)*sigma;
}

void CompactifiedCharacteristicRadialSolution::ingoing_coefficients(const float_type& sigma,
    complex_type& A,complex_type& B,complex_type& C) const {
    coefficients(-1,sigma,A,B,C);
    const complex_type f=horizon_reference_.logarithmic_derivative(sigma);
    B+=2*A*f;
    C=2*constants::I*w_*(1+b_*sigma*sigma)*f-
        (lambda_-float_type(l_)*float_type(l_+1))+2*constants::I*b_*w_*sigma;
}

State CompactifiedCharacteristicRadialSolution::horizon_series(float_type& start) const {
    const float_type d=1-b_,Q=(1+b_)*w_-j_;
    const std::vector<complex_type> A={0,d,3*b_-2,1-3*b_,b_};
    const std::vector<complex_type> P={d-2*constants::I*Q,
        -4+6*b_+4*constants::I*(b_*w_-j_),3-9*b_-2*constants::I*(b_*w_-j_),4*b_};
    const std::vector<complex_type> C={-lambda_-d-2*constants::I*(j_-b_*w_),
        1-3*b_+2*constants::I*(j_-b_*w_),2*b_};
    std::vector<complex_type> h(options_.horizon_order+1,complex_type(0)); h[0]=1;
    for(unsigned n=0;n<options_.horizon_order;++n) {
        complex_type sum=0;
        for(unsigned k=2;k<A.size()&&k<=n+2;++k) {
            unsigned v=n+2-k; sum+=A[k]*float_type(v)*float_type(v==0?0:v-1)*h[v];
        }
        for(unsigned k=1;k<P.size()&&k<=n+1;++k) {
            unsigned v=n+1-k; sum+=P[k]*float_type(v)*h[v];
        }
        for(unsigned k=0;k<C.size()&&k<=n;++k) sum+=C[k]*h[n-k];
        h[n+1]=-sum/(float_type(n+1)*(float_type(n+1)*d-2*constants::I*Q));
    }
    float_type x=float_type("0.1")*minimum(d,1-sigma_max_);
    for(unsigned retry=0;retry<=options_.max_endpoint_retries;++retry) {
        State value=polynomial(h,x); float_type power=1,total=0,tail=0;
        for(unsigned n=0;n<h.size();++n) {
            float_type term=term_size(h[n],n,power); total+=term;
            if(n+4>=h.size()) tail+=term;
            power*=x;
        }
        float_type size=state_size(value,x);
        float_type estimate=tail/maximum(size,float_type("1e-90"));
        const float_type cancellation_indicator=total/maximum(size,float_type("1e-90"));
        float_type roundoff=cancellation_indicator*std::numeric_limits<float_type>::epsilon()*
            float_type(100*(h.size()+1));
        const float_type truncation_estimate=maximum(estimate,roundoff);
        if(truncation_estimate<options_.endpoint_tolerance) {
            start=1-x;value[1]=-value[1];return value;
        }
        x/=2;
    }
    throw std::runtime_error("ENDPOINT_STARTUP_FAILED: horizon Taylor tail");
}

void CompactifiedCharacteristicRadialSolution::static_horizon_polynomials(
    std::vector<complex_type>& A,std::vector<complex_type>& P,
    std::vector<complex_type>& C) const {
    const float_type d=1-b_,Q=(1+b_)*w_-j_;
    const std::vector<complex_type> A0={0,d,3*b_-2,1-3*b_,b_};
    const std::vector<complex_type> P0={d-2*constants::I*Q,
        -4+6*b_+4*constants::I*(b_*w_-j_),3-9*b_-2*constants::I*(b_*w_-j_),4*b_};
    std::vector<complex_type> N(l_+1,complex_type(0));
    const auto& h=horizon_reference_.coefficients();float_type d_power=1;
    for(int k=0;k<=l_;++k) {
        float_type binomial=1;
        for(int n=0;n<=l_-k;++n) {
            N[k+n]+=h[k]/d_power*binomial;
            binomial*=-float_type(l_-k-n)/float_type(n+1);
        }
        d_power*=d;
    }
    const auto sigma_N=product({1,-1},N);
    std::vector<complex_type> Nx(std::max(1,l_),complex_type(0));
    for(int k=1;k<=l_;++k)Nx[k-1]=float_type(k)*N[k];
    auto M=product({-1,1},Nx);M.resize(std::max(M.size(),N.size()),complex_type(0));
    for(std::size_t k=0;k<N.size();++k)M[k]-=float_type(l_+1)*N[k];
    A=product(sigma_N,A0);
    P=product(sigma_N,P0);C=product(sigma_N,
        {-(lambda_-float_type(l_)*float_type(l_+1))+2*constants::I*b_*w_,
         -2*constants::I*b_*w_});
    const auto Pextra=product(A0,M),Cextra=product(
        {2*constants::I*w_*(1+b_),-4*constants::I*w_*b_,2*constants::I*w_*b_},M);
    P.resize(std::max(P.size(),Pextra.size()),complex_type(0));
    C.resize(std::max(C.size(),Cextra.size()),complex_type(0));
    for(std::size_t k=0;k<Pextra.size();++k)P[k]-=2*Pextra[k];
    for(std::size_t k=0;k<Cextra.size();++k)C[k]+=Cextra[k];
}

State CompactifiedCharacteristicRadialSolution::static_horizon_series(float_type& start) const {
    const float_type d=1-b_;
    std::vector<complex_type> A,P,C;
    static_horizon_polynomials(A,P,C);
    std::vector<complex_type> u(options_.horizon_order+1,complex_type(0));u[0]=1;
    for(unsigned n=0;n<options_.horizon_order;++n) {
        complex_type sum=0;
        for(unsigned k=2;k<A.size()&&k<=n+2;++k) {
            unsigned v=n+2-k;sum+=A[k]*float_type(v)*float_type(v==0?0:v-1)*u[v];
        }
        for(unsigned k=1;k<P.size()&&k<=n+1;++k) {
            unsigned v=n+1-k;sum+=P[k]*float_type(v)*u[v];
        }
        for(unsigned k=0;k<C.size()&&k<=n;++k)sum+=C[k]*u[n-k];
        u[n+1]=-sum/(float_type(n+1)*(float_type(n)*A[1]+P[0]));
    }
    float_type x=float_type("0.1")*minimum(d,1-sigma_max_);
    for(unsigned retry=0;retry<=options_.max_endpoint_retries;++retry) {
        State value=polynomial(u,x);float_type power=1,total=0,tail=0;
        for(unsigned n=0;n<u.size();++n) {
            float_type term=term_size(u[n],n,power);total+=term;
            if(n+4>=u.size())tail+=term;power*=x;
        }
        const float_type size=state_size(value,x);
        const float_type cancellation_indicator=total/maximum(size,float_type("1e-90"));
        const float_type truncation_estimate=maximum(tail/maximum(size,float_type("1e-90")),
            cancellation_indicator*std::numeric_limits<float_type>::epsilon()*
            float_type(100*(u.size()+1)));
        if(truncation_estimate<options_.endpoint_tolerance) {
            start=1-x;value[1]=-value[1];return value;
        }
        x/=2;
    }
    throw std::runtime_error("ENDPOINT_STARTUP_FAILED: static horizon Taylor tail");
}

State CompactifiedCharacteristicRadialSolution::infinity_series(float_type& start) const {
    const std::vector<complex_type> A={0,0,1,-(1+b_),b_};
    const std::vector<complex_type> B={-2*constants::I*w_,2,
        -3*(1+b_)-2*constants::I*(b_*w_-j_),4*b_};
    const std::vector<complex_type> C={-lambda_,
        -(1+b_)-2*constants::I*(b_*w_-j_),2*b_};
    std::vector<complex_type> c(options_.infinity_max_order+1,complex_type(0));c[0]=1;
    for(unsigned n=0;n<options_.infinity_max_order;++n) {
        complex_type sum=0;
        for(unsigned k=2;k<A.size()&&k<=n+2;++k) {
            unsigned v=n+2-k;sum+=A[k]*float_type(v)*float_type(v==0?0:v-1)*c[v];
        }
        for(unsigned k=1;k<B.size()&&k<=n+1;++k) {
            unsigned v=n+1-k;sum+=B[k]*float_type(v)*c[v];
        }
        for(unsigned k=0;k<C.size()&&k<=n;++k)sum+=C[k]*c[n-k];
        c[n+1]=-sum/(float_type(n+1)*B[0]);
    }
    float_type s=minimum(float_type("0.01"),minimum(float_type("0.025")*mp::abs(w_),
                                           float_type("0.1")*sigma_min_));
    if(s<=0||s>=sigma_min_)
        throw std::runtime_error("ENDPOINT_STARTUP_FAILED: infinity start must precede source");
    for(unsigned retry=0;retry<=options_.max_endpoint_retries;++retry) {
        std::vector<float_type> terms(c.size()); float_type power=1;
        for(unsigned n=0;n<c.size();++n) {
            terms[n]=term_size(c[n],n,power);power*=s;
        }
        unsigned use=0;
        float_type best=std::numeric_limits<float_type>::max();
        // Use four terms to avoid false convergence at an isolated small coefficient.
        for(unsigned n=8;n+4<c.size();++n) {
            float_type envelope=terms[n+1]+terms[n+2]+terms[n+3]+terms[n+4];
            if(envelope<best) {best=envelope;use=n;}
            if(n>use+8 && envelope>best*100) break;
        }
        std::vector<complex_type> used(c.begin(),c.begin()+use+1);
        State value=polynomial(used,s);
        float_type scale=state_size(value,s),total=0;
        for(unsigned n=0;n<=use;++n)total+=terms[n];
        float_type cancellation=total/maximum(scale,float_type("1e-90"));
        float_type estimate=best/maximum(scale,float_type("1e-90"));
        float_type roundoff=cancellation*std::numeric_limits<float_type>::epsilon()*
            float_type(100*(use+1));
        if(estimate<options_.endpoint_tolerance&&roundoff<options_.endpoint_tolerance) {
            start=s;return value;
        }
        s/=2;
    }
    throw std::runtime_error("ENDPOINT_STARTUP_FAILED: infinity least-term/roundoff estimate");
}

State CompactifiedCharacteristicRadialSolution::static_series(float_type& start) const {
    unsigned order=options_.infinity_max_order;
    std::vector<complex_type> z(order+1,complex_type(0));z[0]=1;
    for(unsigned n=0;n<order;++n) {
        z[n+1]=(float_type(n+l_+1)*((1+b_)*float_type(n+l_+1)-2*constants::I*j_)*z[n]
            -(n?b_*float_type(n+l_)*float_type(n+l_+1)*z[n-1]:complex_type(0)))
            /float_type((n+1)*(n+2*l_+2));
    }
    float_type s=minimum(float_type("0.01"),float_type("0.1")*sigma_min_);
    for(unsigned retry=0;retry<=options_.max_endpoint_retries;++retry) {
        State v=polynomial(z,s);float_type power=1,total=0,tail=0;
        for(unsigned n=0;n<z.size();++n) {
            float_type term=term_size(z[n],n,power);total+=term;
            if(n+4>=z.size())tail+=term;power*=s;
        }
        float_type size=state_size(v,s);
        float_type truncation_estimate=tail/maximum(size,float_type("1e-90"));
        const float_type cancellation_indicator=total/maximum(size,float_type("1e-90"));
        truncation_estimate=maximum(truncation_estimate,
            cancellation_indicator*std::numeric_limits<float_type>::epsilon()*
            float_type(100*(z.size()+1)));
        if(truncation_estimate<options_.endpoint_tolerance) {
            start=s;return v;
        }
        s/=2;
    }
    throw std::runtime_error("ENDPOINT_STARTUP_FAILED: static infinity Taylor tail");
}

void CompactifiedCharacteristicRadialSolution::initialize_endpoints() {
    State in=options_.direct_static_startup?static_horizon_series(in_seed_sigma_):
        horizon_series(in_seed_sigma_);
    if(!options_.direct_static_startup) {
        complex_type S,Ss,Sss;horizon_reference_.evaluate(in_seed_sigma_,S,Ss,Sss);
        in[0]*=in_seed_sigma_/S;
        in[1]=in_seed_sigma_*in[1]/S-(Ss/S-1/in_seed_sigma_)*in[0];
    }
    if(omega_==0)in=State{{complex_type(1),complex_type(0)}};
    // Preserve the additive constant in chi-Omega_H*rstar.
    float_type root=(rp_-rm_)/2;
    float_type chi_minus_oh_rstar=-a_/2-a_/rp_*mp::log(root);
    T_H_=mp::exp(constants::I*m_*chi_minus_oh_rstar)/rp_;
    if(omega_==0&&l_==0&&m_==0) {
        in[0]*=rp_;in[1]*=rp_;T_H_=1;
    }
    in_seed_=in;
    up_seed_=omega_==0?static_series(up_seed_sigma_):infinity_series(up_seed_sigma_);
    T_inf_=1;
    if(omega_==0) {
        integrate_static(-1,in_seed_sigma_,in_seed_);
        integrate_static(1,up_seed_sigma_,up_seed_);
    } else {
        integrate_polynomial_taylor(-1,in_seed_sigma_,in_seed_);
        integrate_polynomial_taylor(1,up_seed_sigma_,up_seed_);
    }
}

void CompactifiedCharacteristicRadialSolution::integrate_static(
    int H,float_type start,const State& initial) {
    auto& dense=H==-1?in_dense_:up_dense_;
    const float_type end=H==-1?sigma_min_:sigma_max_;
    if(start==end) return;
    const float_type direction=H==-1?float_type(-1):float_type(1);
    const float_type t_start=H==-1?float_type(1-start):start;
    const float_type t_end=H==-1?float_type(1-end):end;
    float_type step=minimum(float_type("0.02")*minimum(start,1-start),mp::abs(end-start)/10);
    if(step<=0)throw std::runtime_error("INTEGRATION_FAILED: invalid initial compact step");
    float_type abs_tol=minimum(options_.integration_absolute_tolerance,options_.interpolation_tolerance);
    float_type rel_tol=minimum(options_.integration_relative_tolerance,options_.interpolation_tolerance);
    if(H==1&&options_.outgoing_tolerance>0)
        abs_tol=rel_tol=options_.outgoing_tolerance;
    Stepper stepper(abs_tol,rel_tol,1,1,float_type(0),true);
    stepper.initialize(initial,t_start,step);
    auto rhs=[&](const State& x,State& dx,const float_type& t) {
        const float_type s=H==-1?float_type(1-t):t;
        if(s<=0||s>=1)throw std::runtime_error("INTEGRATION_FAILED: RHS crossed singular endpoint");
        dx[0]=direction*x[1];
        if(H==1) {
            float_type D=(1-s)*(1-b_*s);
            complex_type B=2*(l_+1)-((2*l_+3)*(1+b_)-2*constants::I*j_)*s+
                2*b_*(l_+2)*s*s;
            complex_type C=(l_+1)*(-(l_+1)*(1+b_)+2*constants::I*j_+b_*(l_+2)*s);
            dx[1]=-direction*(B*x[1]+C*x[0])/(s*D);
        } else {
            complex_type A,B,C;
            ingoing_coefficients(s,A,B,C);
            dx[1]=-direction*(B*x[1]+C*x[0])/A;
        }
    };
    std::size_t steps=0;
    while(t_end>stepper.current_time()) {
        if(steps>=options_.max_steps)
            throw std::runtime_error("INTEGRATION_FAILED: compactified step limit");
        float_type remaining=t_end-stepper.current_time();
        if(mp::abs(stepper.current_time_step())>mp::abs(remaining))
            stepper.initialize(stepper.current_state(),stepper.current_time(),remaining);
        stepper.do_step(rhs);++steps;
        const float_type old_sigma=H==-1?float_type(1-stepper.previous_time()):stepper.previous_time();
        const float_type new_sigma=stepper.current_time()==t_end?end:
            (H==-1?float_type(1-stepper.current_time()):stepper.current_time());
        float_type lo=minimum(old_sigma,new_sigma),hi=maximum(old_sigma,new_sigma);
        if(hi>=sigma_min_&&lo<=sigma_max_)dense.emplace_back(lo,hi,stepper);
        if(!mp::isfinite(mp::abs(stepper.current_state()[0]))||
           !mp::isfinite(mp::abs(stepper.current_state()[1])))
            throw std::runtime_error("INTEGRATION_FAILED: nonfinite characteristic state");
    }
    if(H==-1)std::reverse(dense.begin(),dense.end());
}

void CompactifiedCharacteristicRadialSolution::integrate_polynomial_taylor(
    int H,float_type start,const State& initial) {
    auto& dense=H==-1?in_taylor_:up_taylor_;
    const unsigned order=H==-1?32:24;
    const float_type terminal=H==-1?float_type(1-sigma_min_):sigma_max_;
    float_type atol=minimum(options_.integration_absolute_tolerance,options_.interpolation_tolerance);
    float_type rtol=minimum(options_.integration_relative_tolerance,options_.interpolation_tolerance);
    if(H==1&&options_.outgoing_tolerance>0)atol=rtol=options_.outgoing_tolerance;
    std::vector<complex_type> A,B,C;
    if(H==-1) {
        static_horizon_polynomials(A,B,C);
    } else {
        A={0,0,1,-(1+b_),b_};
        B={-2*constants::I*w_,2,-3*(1+b_)-2*constants::I*(b_*w_-j_),4*b_};
        C={-lambda_,-(1+b_)-2*constants::I*(b_*w_-j_),2*b_};
    }
    auto shifted=[](std::vector<complex_type> c,const float_type& center) {
        for(std::size_t n=c.size()-1;n>0;--n)
            for(std::size_t k=n-1;k+1<c.size();++k)c[k]+=center*c[k+1];
        return c;
    };
    const unsigned degree=order+4;
    const float_type eps=std::numeric_limits<float_type>::epsilon()*float_type(100*(degree+1));
    std::size_t steps=0;
    std::vector<float_type> inverse_index(degree+1);
    for(unsigned n=2;n<=degree;++n)inverse_index[n]=float_type(1)/float_type(n*(n-1));
    const float_type infinity_phase_factor=float_type("0.1")*float_type(degree)/mp::abs(w_);
    const float_type horizon_frequency=(1+b_)*w_-j_;
    const float_type horizon_phase_factor=H==-1&&horizon_frequency!=0?
        float_type(float_type("0.1")*float_type(degree)*(1-b_)/
                   (2*mp::abs(horizon_frequency))):float_type(0);
    const float_type radius_factor("0.4");
    const float_type inverse_order=-float_type(1)/order;
    const float_type growth_saturation=mp::pow(float_type("0.4"),order);
    const float_type shrink_saturation=mp::pow(float_type(8),order);
    State state=initial;
    if(H==-1)state[1]=-state[1];
    float_type center=H==-1?float_type(1-start):start;
    float_type step=minimum(float_type("0.2")*center,
        infinity_phase_factor*start*start);
    while(center<terminal) {
        if(steps>=options_.max_steps)
            throw std::runtime_error("INTEGRATION_FAILED: Taylor step limit");
        const auto a=shifted(A,center),b=shifted(B,center),c=shifted(C,center);
        const complex_type inverse_a0=complex_type(1)/a[0];
        std::vector<complex_type> coefficients(degree+1,complex_type(0));
        coefficients[0]=state[0];coefficients[1]=state[1];
        for(unsigned n=0;n+2<=degree;++n) {
            complex_type sum=0;
            for(unsigned k=1;k<a.size()&&k<=n+2;++k) {
                const unsigned v=n+2-k;
                if(v>=2)sum+=a[k]*float_type(v*(v-1))*coefficients[v];
            }
            for(unsigned k=0;k<b.size()&&k<=n+1;++k) {
                const unsigned v=n+1-k;
                if(v)sum+=b[k]*float_type(v)*coefficients[v];
            }
            for(unsigned k=0;k<c.size()&&k<=n;++k)sum+=c[k]*coefficients[n-k];
            coefficients[n+2]=-sum*inverse_a0*inverse_index[n+2];
        }
        std::vector<float_type> magnitudes(degree+1);
        for(unsigned n=0;n<=degree;++n)magnitudes[n]=mp::abs(coefficients[n]);
        step=minimum(step,radius_factor*minimum(center,1-center));
        const float_type sigma=H==-1?float_type(1-center):center;
        step=minimum(step,infinity_phase_factor*sigma*sigma);
        if(horizon_phase_factor!=0)step=minimum(step,horizon_phase_factor*center);
        step=minimum(step,terminal-center);
        State next;
        float_type error=0;
        unsigned retries=0;
        for(;;) {
            if(step<=0||center+step==center)
                throw std::runtime_error("INTEGRATION_FAILED: Taylor step underflow");
            const auto endpoint=polynomial_derivatives(coefficients,step);
            next={{endpoint[0],endpoint[1]}};
            if(!finite_complex(next[0])||!finite_complex(next[1]))
                throw std::runtime_error("INTEGRATION_FAILED: nonfinite Taylor state");
            const float_type value_budget=atol+rtol*maximum(mp::abs(state[0]),mp::abs(next[0]));
            const float_type derivative_budget=atol+rtol*maximum(mp::abs(state[1]),mp::abs(next[1]));
            float_type value_tail=0,derivative_tail=0,value_sum=0,derivative_sum=0,power=1;
            for(unsigned n=0;n<=degree;++n) {
                const float_type derivative_term=n?float_type(magnitudes[n]*float_type(n)*power):float_type(0);
                if(n)power*=step;
                const float_type term=magnitudes[n]*power;
                value_sum+=term;derivative_sum+=derivative_term;
                if(n>order) {
                    value_tail+=term;derivative_tail+=derivative_term;
                }
            }
            error=16*maximum((value_tail+eps*value_sum)/value_budget,
                             (derivative_tail+eps*derivative_sum)/derivative_budget);
            // Check the differentiated interpolant as well as its Taylor tail.
            for(const float_type fraction:{float_type("0.5"),float_type(1)}) {
                const float_type x=fraction*step;
                const auto value=fraction==1?endpoint:polynomial_derivatives(coefficients,x);
                const complex_type aa=polynomial_value(a,x),bb=polynomial_value(b,x),cc=polynomial_value(c,x);
                const float_type acceleration=mp::abs((aa*value[2]+bb*value[1]+cc*value[0])/aa);
                error=maximum(error,16*maximum(acceleration*step*step/value_budget,
                                               acceleration*step/derivative_budget));
            }
            if(error<=1)break;
            if(++retries>40)
                throw std::runtime_error("INTEGRATION_FAILED: Taylor error retry limit");
            step*=error>=shrink_saturation?float_type("0.1"):
                maximum(float_type("0.1"),minimum(float_type("0.8"),
                    float_type("0.8")*mp::pow(error,inverse_order)));
        }
        const float_type end=step==terminal-center?terminal:float_type(center+step);
        if(end>=(H==-1?float_type(1-sigma_max_):sigma_min_))
            dense.push_back({center,end,std::move(coefficients)});
        state=next;center=end;++steps;
        const float_type factor=error<=growth_saturation?float_type(2):
            maximum(float_type("0.5"),minimum(float_type(2),
                float_type("0.8")*mp::pow(error,inverse_order)));
        step*=factor;
    }
}

State CompactifiedCharacteristicRadialSolution::state_at(int H,const float_type& sigma,
                                                         bool characteristic_state) const {
    if(H!=-1&&H!=1)throw std::invalid_argument("characteristic H must be -1 or +1");
    const auto& dense=H==-1?in_dense_:up_dense_;
    const State& seed=H==-1?in_seed_:up_seed_;
    const float_type& seed_sigma=H==-1?in_seed_sigma_:up_seed_sigma_;
    State x;
    if(sigma==seed_sigma)x=seed;
    else if(!(H==-1?in_taylor_:up_taylor_).empty()) {
        const auto& taylor=H==-1?in_taylor_:up_taylor_;
        const float_type coordinate=H==-1?float_type(1-sigma):sigma;
        auto it=std::lower_bound(taylor.begin(),taylor.end(),coordinate,
            [](const TaylorStep& d,const float_type& x){return d.right<x;});
        if(it==taylor.end()||coordinate<it->left||coordinate>it->right)
            throw std::out_of_range("compactified query lies outside Taylor source interval");
        x=polynomial(it->coefficients,coordinate-it->left);
        if(H==-1)x[1]=-x[1];
    } else {
        auto it=std::lower_bound(dense.begin(),dense.end(),sigma,
            [](const DenseStep& d,const float_type& s){return d.right<s;});
        if(it==dense.end()||sigma<it->left||sigma>it->right)
            throw std::out_of_range("compactified query lies outside cached source interval");
        it->snapshot.calc_state(H==-1?float_type(1-sigma):sigma,x);
    }
    if(H==-1&&characteristic_state) {
        complex_type S,Ss,Sss;horizon_reference_.evaluate(sigma,S,Ss,Sss);
        x[1]=S/sigma*x[1]+(Ss/sigma-S/(sigma*sigma))*x[0];x[0]*=S/sigma;
    } else if(H==1&&omega_==0) {
        float_type power=mp::pow(sigma,l_);
        x[1]=power*x[1]+(l_?float_type(l_)*power/sigma*x[0]:complex_type(0));x[0]*=power;
    }
    complex_type scale=H==-1?in_scale_:up_scale_;x[0]*=scale;x[1]*=scale;
    return x;
}

float_type CompactifiedCharacteristicRadialSolution::checked_source_radius(const float_type& r) const {
    if(!mp::isfinite(r))throw std::invalid_argument("radial query requires finite radius");
    // Geodesic turning points may round just outside the configured interval.
    if(r<r_min_) {
        float_type slack=32*std::numeric_limits<float_type>::epsilon()*
            maximum(mp::abs(r_min_),float_type(1));
        if(r_min_-r<=slack)return r_min_;
        throw std::out_of_range("radial query outside configured source interval");
    }
    if(r>r_max_) {
        float_type slack=32*std::numeric_limits<float_type>::epsilon()*
            maximum(mp::abs(r_max_),float_type(1));
        if(r-r_max_<=slack)return r_max_;
        throw std::out_of_range("radial query outside configured source interval");
    }
    return r;
}
void CompactifiedCharacteristicRadialSolution::characteristic(int H,const float_type& r,
    complex_type& y,complex_type& y_sigma) const {
    float_type radius=checked_source_radius(r);
    State x=state_at(H,rp_/radius);y=x[0];y_sigma=x[1];
}
void CompactifiedCharacteristicRadialSolution::reconstruct(int H,const float_type& r,
    const complex_type& y,const complex_type& y_sigma,complex_type& R,complex_type& Rr,
    complex_type& Rrr) const {
    if(H!=-1&&H!=1)throw std::invalid_argument("characteristic H must be -1 or +1");
    if(!mp::isfinite(r)||!finite_complex(y)||!finite_complex(y_sigma)||r<=rp_)
        throw std::invalid_argument("radial reconstruction requires finite exterior data");
    float_type delta=(r-rp_)*(r-rm_),s=rp_/r;
    complex_type G=mp::exp(constants::I*H*(omega_*tortoise(a_,r)-m_*azimuthal_phase(a_,r)))/r;
    complex_type p=constants::I*H*(omega_*(r*r+a_*a_)-a_*m_)/delta-1/r;
    R=G*y;Rr=G*(-s*s/rp_*y_sigma+p*y);
    float_type K=(r*r+a_*a_)*omega_-a_*m_;
    // This reconstruction is not an independent ODE residual.
    Rrr=(-(2*r-2)*Rr-(K*K/delta-lambda_)*R)/delta;
}
void CompactifiedCharacteristicRadialSolution::evaluate_in(const float_type& r,
    complex_type& R,complex_type& Rr,complex_type& Rrr) const {
    float_type radius=checked_source_radius(r);
    const float_type s=rp_/radius,delta=(radius-rp_)*(radius-rm_);
    const State x=state_at(-1,rp_/radius,false);
    const complex_type u=x[0],us=x[1];
    complex_type S,Ss,Sss;horizon_reference_.evaluate(s,S,Ss,Sss);
    const complex_type phase=mp::exp(constants::I*(-omega_*tortoise(a_,radius)+
        m_*azimuthal_phase(a_,radius)))/rp_;
    R=phase*S*u;
    Rr=phase*(-s*s/rp_*(Ss*u+S*us)+
        constants::I*(-omega_*(radius*radius+a_*a_)+a_*m_)/delta*S*u);
    const float_type K=(radius*radius+a_*a_)*omega_-a_*m_;
    Rrr=(-(2*radius-2)*Rr-(K*K/delta-lambda_)*R)/delta;
}
void CompactifiedCharacteristicRadialSolution::evaluate_up(const float_type& r,
    complex_type& R,complex_type& Rr,complex_type& Rrr) const {
    float_type radius=checked_source_radius(r);
    complex_type y,ys;characteristic(1,radius,y,ys);
    reconstruct(1,radius,y,ys,R,Rr,Rrr);
}
complex_type CompactifiedCharacteristicRadialSolution::wronskian(const float_type& r) const {
    float_type radius=checked_source_radius(r);
    complex_type Ri,Rir,Rirr,Ru,Rur,Rurr;
    evaluate_in(radius,Ri,Rir,Rirr);evaluate_up(radius,Ru,Rur,Rurr);
    return (radius-rp_)*(radius-rm_)*(Ri*Rur-Ru*Rir);
}
complex_type CompactifiedCharacteristicRadialSolution::wronskian_characteristic(const float_type& r) const {
    float_type radius=checked_source_radius(r);
    complex_type yi,yis,yu,yus;characteristic(-1,radius,yi,yis);characteristic(1,radius,yu,yus);
    float_type s=rp_/radius,delta=(radius-rp_)*(radius-rm_);
    return (-delta*s*s/rp_*(yi*yus-yu*yis)+
        2*constants::I*(omega_*(radius*radius+a_*a_)-a_*m_)*yi*yu)/(radius*radius);
}
void CompactifiedCharacteristicRadialSolution::rescale_basis(const complex_type& in,
                                                            const complex_type& up) {
    if(!finite_complex(in)||!finite_complex(up)||in==0||up==0)
        throw std::invalid_argument("basis rescaling constants must be finite and nonzero");
    in_scale_*=in;up_scale_*=up;T_H_*=in;T_inf_*=up;
}
