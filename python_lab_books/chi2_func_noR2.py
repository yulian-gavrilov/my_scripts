def Chi2_Jpre_all (Jomegas_pre,Rates,Rates_err):
 	return \
 	(((3./16.*(0*Jomegas_pre[0]+ 1*Jomegas_pre[1]+4*Jomegas_pre[3]))-Rates[0])/Rates_err[0])**2+\
  	(((3./16.*(0*Jomegas_pre[0]+ 3*Jomegas_pre[1]+0*Jomegas_pre[3]))-Rates[2])/Rates_err[2])**2+\
 	(((3./16.*(0*Jomegas_pre[0]+ 1*Jomegas_pre[1]+2*Jomegas_pre[3]))-Rates[3])/Rates_err[3])**2+\
 	(((3./16.*(0*Jomegas_pre[0]+ 1*Jomegas_pre[2]+4*Jomegas_pre[4]))-Rates[4])/Rates_err[4])**2+\
 	(((1./32.*(9*Jomegas_pre[0]+15*Jomegas_pre[2]+6*Jomegas_pre[4]))-Rates[5])/Rates_err[5])**2+\
 	(((3./16.*(0*Jomegas_pre[0]+ 3*Jomegas_pre[2]+0*Jomegas_pre[4]))-Rates[6])/Rates_err[6])**2+\
 	(((3./16.*(0*Jomegas_pre[0]+ 1*Jomegas_pre[2]+2*Jomegas_pre[4]))-Rates[7])/Rates_err[7])**2

