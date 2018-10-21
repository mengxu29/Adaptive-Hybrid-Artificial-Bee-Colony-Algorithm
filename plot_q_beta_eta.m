max_t=10;
[q,beta]=meshgrid(-10:0.2:0.9,0:0.2:10);
eta = (1-q).^(1./beta)*max_t;
mesh(q,beta,eta);