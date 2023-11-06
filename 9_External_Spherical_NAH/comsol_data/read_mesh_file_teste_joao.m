% Rotina para extração de dados de um arquivo .txt -> .mat
% Os arquivos .txt estão organizados de forma sequencial

function [nodes, triangles, data] = read_mesh_file_teste_joao(filename)
% Open the file

fileID = fopen(filename, 'r');

% Read the file line by line
tline = fopen(fileID);
nodes = [];
triangles = [];
data = [];  % 
datas = []; % this will be a column vector that transport '% Data' to [data] in order 
                                              
%     dataIndex = 1;  % Initialize dataIndex
c=1; % contador, will be helpful   

% 

while ischar(tline)
    if c==1 % aplicação de uma condição, as condicionais dentro desse while
        % só irão rodar a primeira vez já que o contador está ==1
        
        if contains(tline, '% Coordinates', 'IgnoreCase', true)
            % Start reading nodes from the next line
            tline = fgetl(fileID);
            while ~contains(tline, '% Elements', 'IgnoreCase', true)
                node = sscanf(tline, '%f');
                nodes = [nodes; node'];
                tline = fgetl(fileID);
            end
        end
        
        if contains(tline, '% Elements (triangles)', 'IgnoreCase', true)
            % Start reading triangles from the next line
            tline = fgetl(fileID);
            while ~contains(tline, 'Data', 'IgnoreCase', true)
                triangle = sscanf(tline, '%d');
                triangles = [triangles; triangle'];
                tline = fgetl(fileID);
            end
        end
        
        if contains(tline, '% Data', 'IgnoreCase', true)
            c=2;
            % Start reading data from the next line
            tline = fgetl(fileID);
            currentData = [];
            while ischar(tline) && ~contains(tline, '% Coordinates', 'IgnoreCase', true)
                datum = sscanf(tline, '%f');
                data = [data; datum'];
                tline = fgetl(fileID);
            end
        end
        tline = fgetl(fileID);

    else % Agora c > 1, a condição imposta por este loop maior não é satisfeita
        % Agora a parte acima é ignorada.
        % A parte abaixo concatena os vetores de cada freq em colunas

        if contains(tline, '% Data', 'IgnoreCase', true)
            % Start reading data from the next line
            tline = fgetl(fileID);
            datas = []; % Zera o vetor de carregamento para cada nova iteração
            while ischar(tline) && ~contains(tline, '% Coordinates', 'IgnoreCase', true)
                datum = sscanf(tline, '%f');
                datas = [datas; datum']; % Vetor de carregamento é carregado
                tline = fgetl(fileID); 
            end
            data(:,c)= datas; % organiza data. Todas as linhas na coluna de valor c
            c = c+1;
        end
        tline = fgetl(fileID);
        
    end
    
end

% Close the file
fclose(fileID);

end
