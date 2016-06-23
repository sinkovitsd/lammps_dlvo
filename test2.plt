dlvo1(d,s,k)=0.5*d*s*s/(2*k*k*1.0/(4*pi*8.987556e6))
dlvo2(d,h)=-0.5*d*h/12.0
dlvo3(d,h,m)=-dlvo2(d,h)*m**8/9.0
oc(x,d,s,k)=dlvo1(d,s,k)*log(1+exp(-k*x))
ov(x,d,h,m)=dlvo2(d,h)/x+dlvo3(d,h,m)/x**9
energy(x,d,s,k,h,m,c,epsr) = (dlvo1(d,s,k)*log(1+exp(-k*x)) - oc(c,d,s,k))/epsr+(dlvo2(d,h)+dlvo3(d,h,m)/x**8)/x-ov(c,d,h,m)
#plot energy(x,2.7,-0.089,1.0/9.6e-2,1.3e-5,0.01,2.5,80)
