z = [2 2]; %partition size
d = [4 4]; %total size
b = d./z;  %number of partitions
N = size(z)(2);
coords_p = [0 0]; %in-partition coordinate contribution
coords_b = [0 0]; %global coordinate contribution


for index = 1:prod(d)
  #begin flat2coords  
  ind_cpy = index;
  part = ceil(index/prod(z));
  div = prod(b(2:N));
  for k = 1:N
    coords_b(k) = (ceil(part/div));
    if (mod(part, div) != 0)
      part = mod(part, div);
    else
      part = div;
    endif
    if (k < N)
      div = div/b(k+1);
    endif
  endfor
  #end flat2coords pt 1
  
  
  #localize the index (do in ravel proper!)
  if (mod(ind_cpy,prod(z)) != 0)
    ind_cpy = mod(ind_cpy, prod(z));
  else
    ind_cpy = prod(z);
  endif
  #begin flat2coords call 2)
  div = prod(z(2:N));  
  for k = 1:N
    ##      if (mod(ceil(index/div),z(k)) != 0)
    ##      coords_p(k) = mod(ceil(index/div), z(k));
    ##    else
    ##      coords_p(k) = z(k);
    ##    endif
    coords_p(k) = ceil(ind_cpy/div);
    if (mod(ind_cpy, div) != 0)
      ind_cpy = mod(ind_cpy, div);
    else
      ind_cpy = div;
    endif
    
    if (k < N)
      div = div/z(k + 1);
    endif
    #end flat2coords call 2
  endfor
  
  coords_p;
  coords_b;
  coords = coords_p + ((coords_b)-1).*z;
  
  disp(coords)
  #end flat2coords
  
  #begin coords2flat
  flat = coords(N);
  div = prod(d);
  for k = 1:N-1
    div = div/d(k);
    flat = flat + (coords(k) - 1)*div;
  endfor
  #end coords2flat
  
  #disp(ind_cpy)
  disp(flat)
  endfor