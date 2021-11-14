ZZ11 = B11 + skel_struct.U11*skel_struct.V11;
ZZ12 = skel_struct.U12*skel_struct.V12;
ZZ21 = skel_struct.U21*skel_struct.V21;
ZZ22 = B22 + skel_struct.U22*skel_struct.V22;

ZZ = [ZZ11,ZZ12;ZZ21,ZZ22];