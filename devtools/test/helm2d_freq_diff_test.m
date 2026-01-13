%%%%%%%%%% TEST FOR THE exponential %%%%%%%%%%%

clear;
% zk = 1:0.1:100;
% err = zeros(1,length(zk));
% h0 = 1e-3;
% for i = 1:length(zk)
%     h = h0/zk(i);
%     xp = exp(pi + h);
%     xm = exp(pi - h);
%     derx = exp(pi);
%     err(i) = abs((1/h)*(xp - xm) - 2*derx);
% end
% loglog(zk, err, 'k.');
% hold on; loglog(zk, 8*(h0./zk).^2, 'b.');
% return

%%%%%%%%%% TEST FOR THE freq_diff %%%%%%%%%%%

chnkr = get_chunker(1);
zk = 1:0.1:30;
err = zeros(1,length(zk));
err1 = zeros(1,length(zk));
h0 = 1e-4;
zz = 1.2;
for i = 1:length(zk)
     h = h0/zk(i);
     Dp = kernel('helm', 'd', zz + h);  
     Ap = chunkermat(chnkr, Dp);
     Dm = kernel('helm', 'd', zz - h);
     Am = chunkermat(chnkr, Dm);
     Freqk = zz*kernel('helm', 'freq_diff', zz);
     Ak = chunkermat(chnkr, Freqk); 
     err1(i) = norm(Ap - Am);
     % err(i) = norm((1/h)*(Ap - Am) - 2*Ak);
     err(i) = norm(Ap - Am - 2*h*Ak);
end
loglog(zk, err, 'k.');
hold on; loglog(zk, (h0./zk).^3, 'b.');
%hold on; loglog(zk, (h0./zk), 'r.');



% zk = 0:0.01:1;
% Dk = kernel('helm', 'd', 3);  
% Freqk = kernel('helm', 'freq_diff', 3);
% Ak = chunkermat(chnkr, Dk);  
% Fk = chunkermat(chnkr, Freqk);
% for i = 1:length(zk)
%     Dp = kernel('helm', 'd', 3 + zk(i));  
%     Ap = chunkermat(chnkr, Dk); 
%     Dm = kernel('helm', 'd', 3 - zk(i));
%     Am = chunkermat(chnkr, Dk);
%     Freq = kernel('helm', 'freq_diff', 3 + zk(i));
%     F = chunkermat(chnkr, Freqk);
% end
% %%
% 
%     det1 = det(A);
% end
