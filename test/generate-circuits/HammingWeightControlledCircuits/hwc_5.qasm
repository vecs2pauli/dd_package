// This circuit is a randomly generated Hamming Weight-controlled circuit. The control and target registers have 5 qubits each

OPENQASM 2.0;
include "qelib1.inc";
qreg q[17];

// Step 1: put the first register in uniform superposition, by applying hadamard to all qubits
h q[16];
h q[15];
h q[14];
h q[13];
h q[12];
// Stage 1 of Hamming control gates
// Stage 1: compute part  (5 qubits)
x q[16];
x q[15];
ccx q[16], q[15], q[11];
x q[16];
x q[15];
x q[16];
ccx q[16], q[15], q[10];
x q[16];
x q[15];
ccx q[16], q[15], q[10];
x q[15];
ccx q[16], q[15], q[9];
// Stage 2 of Hamming control gates
// Stage 2: compute part  (5 qubits)
x q[14];
ccx q[14], q[11], q[5];
x q[14];
ccx q[14], q[11], q[4];
x q[14];
ccx q[14], q[10], q[4];
x q[14];
ccx q[14], q[10], q[3];
x q[14];
ccx q[14], q[9], q[3];
x q[14];
ccx q[14], q[9], q[2];
// Stage 2: uncompute part  (n=5 qubits)
x q[16];
x q[15];
ccx q[16], q[15], q[11];
x q[16];
x q[15];
x q[16];
ccx q[16], q[15], q[10];
x q[16];
x q[15];
ccx q[16], q[15], q[10];
x q[15];
ccx q[16], q[15], q[9];
// Stage 2: swap part  (5 qubits)
//swap q[11],q[5];
cx q[11],q[5];
cx q[5],q[11];
cx q[11],q[5];
//swap q[10],q[4];
cx q[10],q[4];
cx q[4],q[10];
cx q[10],q[4];
//swap q[9],q[3];
cx q[9],q[3];
cx q[3],q[9];
cx q[9],q[3];
//swap q[8],q[2];
cx q[8],q[2];
cx q[2],q[8];
cx q[8],q[2];
// Stage 3 of Hamming control gates
// Stage 3: compute part  (5 qubits)
x q[13];
ccx q[13], q[11], q[5];
x q[13];
ccx q[13], q[11], q[4];
x q[13];
ccx q[13], q[10], q[4];
x q[13];
ccx q[13], q[10], q[3];
x q[13];
ccx q[13], q[9], q[3];
x q[13];
ccx q[13], q[9], q[2];
x q[13];
ccx q[13], q[8], q[2];
x q[13];
ccx q[13], q[8], q[1];
// Stage 3: uncompute part  (5 qubits)
Cx q[5], q[11];
x q[13];
ccx q[13], q[4], q[10];
x q[13];
ccx q[13], q[4], q[9];
x q[13];
ccx q[13], q[3], q[9];
x q[13];
ccx q[13], q[3], q[8];
x q[13];
ccx q[13], q[2], q[8];
x q[13];
ccx q[13], q[2], q[7];
ccx q[13], q[1], q[8];
// Stage 3: swap part  (5 qubits)
//swap q[11],q[5];
cx q[11],q[5];
cx q[5],q[11];
cx q[11],q[5];
//swap q[10],q[4];
cx q[10],q[4];
cx q[4],q[10];
cx q[10],q[4];
//swap q[9],q[3];
cx q[9],q[3];
cx q[3],q[9];
cx q[9],q[3];
//swap q[8],q[2];
cx q[8],q[2];
cx q[2],q[8];
cx q[8],q[2];
//swap q[7],q[1];
cx q[7],q[1];
cx q[1],q[7];
cx q[7],q[1];
// Stage 4 of Hamming control gates
// Stage 4: compute part  (5 qubits)
x q[12];
ccx q[12], q[11], q[5];
x q[12];
ccx q[12], q[11], q[4];
x q[12];
ccx q[12], q[10], q[4];
x q[12];
ccx q[12], q[10], q[3];
x q[12];
ccx q[12], q[9], q[3];
x q[12];
ccx q[12], q[9], q[2];
x q[12];
ccx q[12], q[8], q[2];
x q[12];
ccx q[12], q[8], q[1];
x q[12];
ccx q[12], q[7], q[1];
x q[12];
ccx q[12], q[7], q[0];
// Stage 4: uncompute part  (5 qubits)
Cx q[5], q[11];
x q[12];
ccx q[12], q[4], q[10];
x q[12];
ccx q[12], q[4], q[9];
x q[12];
ccx q[12], q[3], q[9];
x q[12];
ccx q[12], q[3], q[8];
x q[12];
ccx q[12], q[2], q[8];
x q[12];
ccx q[12], q[2], q[7];
x q[12];
ccx q[12], q[1], q[7];
x q[12];
ccx q[12], q[1], q[6];
ccx q[12], q[0], q[7];
// Stage 4: swap part  (5 qubits)
//swap q[11],q[5];
cx q[11],q[5];
cx q[5],q[11];
cx q[11],q[5];
//swap q[10],q[4];
cx q[10],q[4];
cx q[4],q[10];
cx q[10],q[4];
//swap q[9],q[3];
cx q[9],q[3];
cx q[3],q[9];
cx q[9],q[3];
//swap q[8],q[2];
cx q[8],q[2];
cx q[2],q[8];
cx q[8],q[2];
//swap q[7],q[1];
cx q[7],q[1];
cx q[1],q[7];
cx q[7],q[1];
//swap q[6],q[0];
cx q[6],q[0];
cx q[0],q[6];
cx q[6],q[0];
// The Hamming-weight register has now been set up.
// Next, we add 5 Hadamard gates and 95 random controlled Clifford gates
h q[5];
h q[4];
h q[3];
h q[2];
h q[1];
// Next we add the 95 random controlled Clifford gates
ccz q[7],q[2],q[4];
ccx q[10],q[1],q[4];
cx q[6],q[1];
cy q[9],q[2];
ccy q[9],q[2],q[3];
ccy q[8],q[4],q[1];
ccy q[6],q[2],q[5];
ccz q[11],q[1],q[2];
ccz q[9],q[1],q[3];
cx q[8],q[5];
ccx q[10],q[1],q[3];
ccz q[6],q[2],q[1];
ccy q[8],q[4],q[3];
ch q[10],q[3];
cx q[8],q[2];
ccx q[7],q[2],q[5];
ccz q[9],q[5],q[4];
cs q[8],q[1];
cs q[10],q[5];
ccz q[7],q[3],q[4];
cx q[7],q[5];
cs q[8],q[3];
ch q[11],q[5];
ccx q[6],q[5],q[1];
ccx q[8],q[4],q[1];
ccy q[11],q[2],q[5];
ccx q[8],q[1],q[2];
cs q[8],q[4];
ccz q[10],q[1],q[4];
ch q[10],q[4];
ccz q[11],q[1],q[2];
ccx q[10],q[5],q[1];
ch q[11],q[5];
ch q[6],q[2];
cz q[10],q[1];
ccy q[9],q[5],q[3];
ccy q[8],q[3],q[2];
ccy q[6],q[3],q[1];
ccx q[6],q[5],q[2];
cz q[8],q[5];
ccy q[8],q[3],q[2];
ccx q[6],q[3],q[4];
cx q[8],q[5];
ccx q[6],q[2],q[3];
cy q[8],q[3];
ccy q[8],q[4],q[5];
ch q[7],q[1];
cy q[9],q[4];
ccx q[8],q[1],q[2];
ccz q[11],q[5],q[4];
cz q[9],q[1];
ccy q[6],q[5],q[3];
ccz q[11],q[5],q[3];
cx q[7],q[1];
cz q[11],q[5];
cs q[6],q[4];
ch q[7],q[4];
cs q[8],q[4];
ccx q[9],q[5],q[2];
ch q[10],q[1];
cx q[10],q[5];
cz q[6],q[2];
ccz q[6],q[1],q[3];
ccz q[9],q[2],q[5];
cz q[8],q[3];
ccx q[7],q[4],q[5];
cy q[8],q[5];
ccy q[10],q[5],q[2];
ccz q[6],q[5],q[1];
cz q[11],q[4];
cx q[7],q[2];
ccz q[9],q[2],q[1];
cs q[7],q[5];
ch q[11],q[4];
ccz q[9],q[2],q[1];
ccx q[6],q[1],q[3];
ccx q[7],q[2],q[5];
ccy q[7],q[1],q[5];
ccy q[11],q[1],q[3];
ch q[7],q[2];
cy q[10],q[4];
ch q[10],q[5];
ccy q[7],q[2],q[1];
cx q[9],q[4];
cy q[6],q[5];
cz q[11],q[4];
ccz q[11],q[5],q[3];
ccy q[7],q[2],q[3];
ch q[9],q[2];
cz q[6],q[4];
cz q[8],q[2];
ccx q[9],q[3],q[1];
ch q[7],q[5];
cy q[8],q[1];
ccx q[9],q[3],q[1];
