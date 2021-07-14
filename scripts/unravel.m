z = [2 3]; %partition size
d = [4 6]; %total size
b = d./z; % number of partitions
N = size(z)(2);
coords_p = [0 0]; %in-partition coordinate contribution
coords_b = [0 0]; %global coordinate contribution


for index = 1:prod(d)
  #begin flat2coords  
  ind_cpy = index;
  div = prod(b)/b(1);
  part = ceil(index/prod(z));
  
  div = prod(b)/b(N);
    for k = 1:N
    if (mod(ceil(part/div),b(k)) != 0)
      coords_b(k) = mod(ceil(part/div), b(k));
    else
      coords_b(k) = b(k);
    endif
    
    if (k > 1)
      div = div/b(k);
    endif
    
  endfor
  
  div = prod(z)/z(N);
  
  for k = 1:N
    if (mod(ceil(index/div),z(k)) != 0)
      coords_p(k) = mod(ceil(index/div), z(k));
    else
      coords_p(k) = z(k);
    endif
    
    if (k > 1)
      div = div/z(k);
    endif
    
  endfor
 
  coords_p;=
  coords_b;
  coords = coords_p + ((coords_b)-1).*z
  #end flat2coords
  
  #begin coords2flat
  flat = coords(1);
  for k = 2:N
    flat = flat + (coords(k) - 1)*prod(d(1:k-1));
  endfor
  #end coords2flat
  
  #disp(ind_cpy)
  disp(flat)
endfor