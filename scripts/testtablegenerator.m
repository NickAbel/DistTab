table_dims = [6, 6]; %[CVAR, CMEAN]
CVAR = linspace(0,1,table_dims(1));
CMEAN = linspace(0,1,table_dims(2));
block_dims = [3, 3];
total_blocks = table_dims./block_dims;

file_id = fopen('test_block_table.dat','w');

fprintf(file_id, '%d %d\n', table_dims(1), table_dims(2))
fprintf(file_id, '%d %d\n', block_dims(1), block_dims(2))
m = 0;

block_counter = 1;
for i = 1:total_blocks(2)
  for j = 1:total_blocks(1)
    for k = 1:block_dims(2)
      for l = 1:block_dims(1)
        m = m + 1;
       fprintf(file_id, '%.7f %.7f %.7f\n',CVAR(l + (j-1)*block_dims(1)), CMEAN(k + (i-1)*block_dims(2)), block_counter + m/1000.0)
      endfor
    endfor
    block_counter = block_counter + 1;
    m = 0;
  endfor
endfor

fclose(file_id);