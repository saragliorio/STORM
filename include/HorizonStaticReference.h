#pragma once

#include "config.h"
#include <stdexcept>
#include <vector>

class HorizonStaticReference {
public:
    explicit HorizonStaticReference(int ell, float_type b=0, float_type j=0) : d_(1-b) {
        if(ell<0||!mp::isfinite(b)||!mp::isfinite(j)||b<0||b>=1)
            throw std::invalid_argument("invalid static horizon reference");
        h_.resize(static_cast<std::size_t>(ell)+1);h_[0]=1;
        const complex_type c=1+2*constants::I*j/d_;
        for(int k=0;k<ell;++k)
            h_[k+1]=h_[k]*float_type(ell+k+1)*float_type(ell-k)/
                (float_type(k+1)*(c+k));
    }
    const std::vector<complex_type>& coefficients() const { return h_; }
    void evaluate(const float_type& sigma, complex_type& S, complex_type& S_sigma,
        complex_type& S_sigmasigma) const {
        if(!mp::isfinite(sigma)||sigma<=0||sigma>1)
            throw std::invalid_argument("static horizon reference requires 0<sigma<=1");
        const float_type q=(1-sigma)/(d_*sigma),qs=-1/(d_*sigma*sigma);
        S=h_.back();S_sigma=0;S_sigmasigma=0;
        for(std::size_t k=h_.size()-1;k>0;--k) {
            S_sigmasigma=S_sigmasigma*q+2*S_sigma;
            S_sigma=S_sigma*q+S;S=S*q+h_[k-1];
        }
        S_sigmasigma=S_sigmasigma*qs*qs-2*qs/sigma*S_sigma;
        S_sigma*=qs;
        if(S==0||!mp::isfinite(S.real())||!mp::isfinite(S.imag())||
           !mp::isfinite(S_sigma.real())||!mp::isfinite(S_sigma.imag())||
           !mp::isfinite(S_sigmasigma.real())||!mp::isfinite(S_sigmasigma.imag()))
            throw std::runtime_error("STATIC_REFERENCE_FAILED: zero or nonfinite polynomial");
    }
    complex_type logarithmic_derivative(const float_type& sigma) const {
        complex_type S,Ss,Sss;evaluate(sigma,S,Ss,Sss);return Ss/S-1/sigma;
    }
private:
    float_type d_;
    std::vector<complex_type> h_;
};
