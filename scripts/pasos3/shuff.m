

function vect_shuff=shuff(vect)
% shuffle the column vector vect producing a randomly permuted version built 
%to control against the null-hiphothesis, of being by chance
n=size(vect,1);
[~,idx]=sort(rand(n,1));
vect_shuff=vect(idx);

end