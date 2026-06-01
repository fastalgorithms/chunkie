stokes_dtracTest0();


function stokes_dtracTest0()
srcinfo = [];
srcinfo.r = rand(2,1);
srcinfo.n = rand(2,1);

targinfo = [];
targinfo.r = rand(2,1);
targinfo.n = rand(2,1);

strengths = rand(2,1);

mu = 1.1;

kernt = kernel('stok', 'dtrac', mu);
kerng = kernel('stok', 'dgrad', mu);
kernp = kernel('stok', 'dpres', mu);

Kt = kernt.eval(srcinfo, targinfo)*strengths;
Kg = kerng.eval(srcinfo, targinfo)*strengths;
Kp = kernp.eval(srcinfo, targinfo)*strengths;

du = reshape(Kg, [2,2,1]);
dut = permute(du, [2,1,3]);
eu = du + dut;

euxx = squeeze(eu(1,1,:));
euxy = squeeze(eu(1,2,:));
euyy = squeeze(eu(2,2,:));
f = zeros(2,1);
p = Kp.';
ntx = targinfo.n(1,:).';
nty = targinfo.n(2,:).';
f(1:2:end) = -p.*ntx + (euxx.*ntx + euxy.*nty)*mu;
f(2:2:end) = -p.*nty + (euxy.*ntx + euyy.*nty)*mu;

norm(f - Kt)


end


