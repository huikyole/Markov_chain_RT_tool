function [P11_nr]=HG_Phase(g,costheta)
P11_nr=(1-g^2)./[(1+g^2-2*g*costheta).^1.5];
end