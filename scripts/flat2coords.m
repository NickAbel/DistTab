z = [2 2]; %partition size
d = [4 4]; %total size
b = d./z; % number of partitions
N = size(z)(2);
coords_p = [0 0]; %in-partition coordinate contribution
coords_b = [0 0]; %global coordinate contribution


for index = 1:prod(d)
  ind_cpy = index;
  div = prod(b)/b(N);
  part = ceil(index/prod(z));
  
  div = prod(b)/b(N);
    for k = 1:N
    if (mod(ceil(part/div),b(N - k + 1)) != 0)
      coords_b(N - k + 1) = mod(ceil(part/div), b(N - k + 1));
    else
      coords_b(N - k + 1) = b(N - k + 1);
    endif
    
    if (mod(part,div) != 0)
      part = mod(part, div);
    else
      part = div;
    endif
    
    if (k < N)
      div = div/b(N - k);
    endif
    
  endfor
  
  div = prod(z)/z(N);
  
  if (mod(ceil(index/div),z(N - k + 1)) != 0)
    coords_p(N - k + 1) = mod(ceil(index/div), z(N - k + 1));
  else
    coords_p(N - k + 1) = z(N - k + 1);
  endif
  
  for k = 1:N
    if (mod(ceil(index/div),z(N - k + 1)) != 0)
      coords_p(N - k + 1) = mod(ceil(index/div), z(N - k + 1));
    else
      coords_p(N - k + 1) = z(N - k + 1);
    endif
    
    if (mod(index,div) != 0)
      index = mod(index, div);
    else
      index = div;
    endif
    
    if (k < N)
      div = div/z(N - k);
    endif
    
  endfor
  coords_p;
  coords_b;
  coords = coords_p + ((coords_b)-1).*z;
  flat = coords(1);
  for k = 2:N
    flat = flat + (coords(k) - 1)*prod(d(1:k-1));
  endfor
  #disp(ind_cpy)
  disp(flat)
endfor