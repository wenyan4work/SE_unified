function [ f ] = FullErrFunc( varkc, L, xi )

f=...
(16.*pi.^2.*((3.*varkc.^2 + 2.*(3 + (L.^2.*varkc.^2)./2.).*xi.^2)./ ...
exp(varkc.^2./(4..*xi.^2)) +  ...
((0+1i.*1).*L.^4.*sqrt(pi).*varkc.*xi.^5.* ...
((0+1i.*2) - (0+1i.*1).*erfz(varkc./(2..*xi) - (0+1i.*1).*L.*xi) -  ...
(0+1i.*1).*erfz(varkc./(2..*xi) + (0+1i.*1).*L.*xi)))./exp(L.^2.*xi.^2)))./ ...
(varkc.*xi.^2);

end

