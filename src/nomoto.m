function nomoto(T1,T2,T3,K)
% NOMOTO(T1,T2,T3,K) generates the Bode plots for
%
%
% Author:Thor I. Fossen

T = T1+T2-T3;
d1 = [T 1 0];
n1 = K;
d2 = [T1*T2 T1+T2 1 0];
n2 = K*[T3 1];
[mag1,phase1,w] = bode(n1,d1);
[mag2,phase2] = bode(n2,d2,w);
% shift ship phase with 360 deg for course unstable ship
if K < 0
    phase1 = phase1-360;
    phase2 = phase2-360;
end
clf,subplot(211),semilogx(w,20*log10(mag1))
grid
xlabel('Frequency [rad/s]'),title('Gain [dB]')
hold on,semilogx(w,20*log10(mag2),'–'),hold off
subplot(212),semilogx(w,phase1),grid
xlabel('Frequency [rad/s]'),title('Phase [deg]')
hold on,semilogx(w,phase2,'–'),hold off
end
