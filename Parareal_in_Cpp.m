void parareal (float T, float u0 [nsta], int N, int K, float U[10] [N+1] [nsta], int MG, int MF)
{float dT = TN; float TT [N]; int i;
for(i = 0; i < N+1; i++) {
TT[i] = i * dT; }
float Go [N+1] [nsta]; float Gn [N+1] [nsta]; float Fn [N+1] [nsta]; 
for (i = 0; i < nsta; i++) {
U[0][0][i]= u0 [i]; }
int n;
for ( n = 0; n < N; n++) {
G(TT [n], TT [n + 1], U[0] [n], MG, Go [n + 1]);
for (i=0; i<nsta; i++) {
U[0] [n 1] [i] = Go [n + 1] [i]; } }
int k; int m; int A = 1000; int B = 1;
for (k 0; k < K; k++) {
for ( n = 0; n < A; n++) {
#pragma omp parallel private (m)
#pragma omp for
for (m = 0; mB; m++) {
#OpenMP Code
F(TT[m + B*n ], TT [m + B*n + 1], U[k] [m+B*n], MF, Fn [m + B*n + 1]); }}}
for (i = 0; i<nsta; i++) {
Uk 1][0][i] = u0 [i]; }
for ( n = 0; n < N; n++) {
G(TT [n], TT [n + 1], U[k+1] [n], MG, Gn [n + 1]);
for (i=0; i<nsta; i++) {
U[k + 1] [n + 1] [i] = Fn [n+1] [i] + Gn [n+1] [i] - Go [n+1] [i]; } } 
for ( n = 0; n < N; n++) 
    { for (i=0; i<nsta; i++) {
Go [n] [i] = Gn [n] [i]; 
}}}}
