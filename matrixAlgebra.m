
printf("Example matrix A\n");
A = [ 24790661.71, -296708.9463, 52630.21512; -296708.9463,  4637785.923, -918282.2373; 52630.21512, -918282.2373, 4728995.838 ]
printf("Get eigenvalues and eigenvectors\n");
[V, lambda] = eig(A);
printf("Eigenvalues:\n ");
lambda
printf("Eigenvectors:\n ");
V
printf("Show that eigenvectors are orthonormal V' = inv(V)\n");
printf("V'\nc");
V'
printf("invV\n");
inv(V)
printf("V'*V\n");
V'*V
% reconstruct A given Q and V
V*lambda*inv(V)  % this is correct


