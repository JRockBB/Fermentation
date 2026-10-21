functions {
	real[] 	cbf(real t,
				real[] x,
				real[] theta,
				real[] y_r,
				int[] y_i) {
		real dxdt[8];
		real v1;
		real v2;
		real v3;
		real v4;
		real v5;
		real v6;
		real v7;
		real v8;
		real v9;
		real v10;
		real v11;
		
		v1 = theta[12]*x[1]*x[6]/(theta[17]+x[1]);
		v2 = theta[13]*x[2]*x[6]/(theta[18]+x[2]);
		v3 = theta[14]*x[1]*x[7]/(theta[19]+x[1]);
		v4 = theta[15]*x[3]*x[8]/(theta[20]+x[3]);
		v5 = theta[16]*x[4]*x[8]/(theta[21]*x[8]+x[4]);
		v6 = theta[22]*x[6]*x[3];
		v7 = theta[23]*x[7]*x[4];
		v8 = theta[24]*x[8]*x[5]^2;
		v9 = theta[25]*x[2]*x[7]/(theta[26]+x[2]);
		v10 = theta[34]*x[4]*x[6]/(theta[35]+x[4]);
		v11 = theta[38]*x[5]*x[8]/(theta[39]+x[5]);
		
		dxdt[1] = - theta[1]*v1 - theta[2]*v3;
		dxdt[2] = - theta[3]*v2 - theta[27]*v9;
		dxdt[3] = theta[4]*v1 + theta[5]*v2 + theta[28]*v3 + theta[29]*v9 - theta[6]*v4 + theta[36]*v10;
		dxdt[4] = theta[7]*v3 + theta[30]*v9 - theta[8]*v5 - theta[37]*v10;
		dxdt[5] = theta[9]*v3 + theta[31]*v9 + theta[10]*v4 + theta[11]*v5 + theta[32]*v1 + theta[33]*v2 - theta[40]*v11;
		dxdt[6] = v1 + v2 - v6 + v10;
		dxdt[7] = v3 + v9 - v7;
		dxdt[8] = v4 + v5 - v8 + v11;
		
		return dxdt;
	}
}

data {
	int<lower=1> T;
	real<lower=0> x[T,8];
	real t0;
	real ts[T];
	real x0[8];
	real scl[8];
}

transformed data {

	real y_r[0];
	int y_i[0];
	real<lower=0> x0_1[8];
	real<lower=0> xn[T,8];
	for (t in 1:T)
		for (n in 1:8)
			xn[t,n]=x[t,n]/scl[n];
  
  for (n in 1:8)
    x0_1[n] = x0[n]/scl[n]; 
}
parameters {
	real <lower=0> sigma;
	real <lower=0> mu1;
	real <lower=0> mu2;
	real <lower=0> mu3;
	real <lower=0> mu4;
	real <lower=0> mu5;
	real <lower=0> mu6;
	real <lower=0> ks1;
	real <lower=0> ks2;
	real <lower=0> ks3;
	real <lower=0> ks4;
	real <lower=0> ks5;
	real <lower=0> ks6;
	real <lower=0> k1;
	real <lower=0> k2;
	real <lower=0> k3;
	real <lower=0> yc1;
	real <lower=0> yc2;
	real <lower=0> yc3;
	real <lower=0> yc4;
	real <lower=0> yc5;
	real <lower=0> yc6;
	real <lower=0> yc7;
	real <lower=0> yc8; 
	real <lower=0> yc9;
	real <lower=0> yc10;
	real <lower=0> yc11;
	real <lower=0> yc12;
	real <lower=0> yc13;
	real <lower=0> yc14;
	real <lower=0> yc15;
	real <lower=0> yc16;
	real <lower=0> yc17;
	real <lower=0> yc18;
	real <lower=0> mu7;
	real <lower=0> ks7;
	real <lower=0> yc19;
	real <lower=0> yc20;
	real <lower=0> mu8;
	real <lower=0> ks8;
	real <lower=0> yc21;
}

transformed parameters {
	real x_hat[T,8];
	{
		real theta[40];
		theta[1] = yc1;
		theta[2] = yc2;
		theta[3] = yc3;
		theta[4] = yc4;
		theta[5] = yc5;
		theta[6] = yc6;
		theta[7] = yc7;
		theta[8] = yc8;
		theta[9] = yc9;
		theta[10] = yc10;
		theta[11] = yc11;
		theta[12] = mu1;
		theta[13] = mu2;
		theta[14] = mu3;
		theta[15] = mu4;
		theta[16] = mu5;
		theta[17] = ks1;
		theta[18] = ks2;
		theta[19] = ks3;
		theta[20] = ks4;
		theta[21] = ks5;
		theta[22] = k1;
		theta[23] = k2;
		theta[24] = k3;
		theta[25] = mu6;
		theta[26] = ks6;
		theta[27] = yc12;
		theta[28] = yc13;
		theta[29] = yc14;
		theta[30] = yc15;
		theta[31] = yc16;
		theta[32] = yc17;
		theta[33] = yc18;
		theta[34] = mu7;
		theta[35] = ks7;
		theta[36] = yc19;
		theta[37] = yc20;
		theta[38] = mu8;
		theta[39] = ks8;
		theta[40] = yc21;

		x_hat = integrate_ode_rk45(cbf, x0_1, t0, ts, theta, y_r, y_i,1.0E-6, 1.0E-6, 1.0E4);
	}
}
model{
	mu1~normal(0.5,0.3);
	mu2~normal(0.5,0.3);
	mu3~normal(0.5,0.3);
	mu4~normal(0.5,0.3);
	mu5~normal(0.5,0.3);
	mu6~normal(0.5,0.3);
	ks1~normal(0.5,0.3);
	ks2~normal(0.5,0.3);
	ks3~normal(0.5,0.3);
	ks4~normal(0.5,0.3);
	ks5~normal(0.5,0.3);
	ks6~normal(0.5,0.3);
	k1~normal(0.5,0.3);
	k2~normal(0.5,0.3);
	k3~normal(0.5,0.3);
	yc1~normal(0.5,0.3);
	yc2~normal(0.5,0.3);
	yc3~normal(0.5,0.3);
	yc4~normal(0.5,0.3);
	yc5~normal(0.5,0.3);
	yc6~normal(0.5,0.3);
	yc7~normal(0.5,0.3);
	yc8~normal(0.5,0.3);
	yc9~normal(0.5,0.3);
	yc10~normal(0.5,0.3);
	yc11~normal(0.5,0.3);
	yc12~normal(0.5,0.3);
	yc13~normal(0.5,0.3);
	yc14~normal(0.5,0.3);
	yc15~normal(0.5,0.3);
	yc16~normal(0.5,0.3);
	yc17~normal(0.5,0.3);
	yc18~normal(0.5,0.3);
	mu7~normal(0.5,0.3);
	ks7~normal(0.5,0.3);
	yc19~normal(0.5,0.3);
	yc20~normal(0.5,0.3);
	mu8~normal(0.5,0.3);
	ks8~normal(0.5,0.3);
	yc21~normal(0.5,0.3);

	sigma~cauchy(0,1);
	for (t in 1:T)
		xn[t]~normal(x_hat[t],sigma);
}

generated quantities{
  real log_lik[T,8];
  for (t in 1:T)
    for (n in 1:8)
      log_lik[t,n] = normal_lpdf(xn[t,n]|x_hat[t,n],sigma);
}

