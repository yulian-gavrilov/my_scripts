def Chi2_Jpre_all (Jomegas_pre,Rates,Rates_err):
 	return \
       	(((3./16.*(0*Jomegas_pre[0]+ 3*Jomegas_pre[2]+0*Jomegas_pre[4]))-Rates[6])/Rates_err[6])**2+\
 	(((3./16.*(0*Jomegas_pre[0]+ 1*Jomegas_pre[2]+2*Jomegas_pre[4]))-Rates[7])/Rates_err[7])**2

