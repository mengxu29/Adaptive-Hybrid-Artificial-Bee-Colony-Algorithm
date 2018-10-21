%% Intervalos de confianca: bootstrap parametrico
B = 999 ;                                                                    % Numero de amostras bootstrap.
n = size(t(:,1),1) ;                                                            % Tamanho da amostra.
boot = zeros(B+1,4) ;                                                       % Matriz de parametros [eta beta q] com dimensao B x 3.
                                             % Primeiro elemento da matriz e a solucao da amostra original.

eta = 4455.2041378061148862 ;%mPart.p(1,1) ; % 
beta = 0.7725982218311025 ; %mPart.p(1,2) ; % 
q = -2.1910149665244103 ;%mPart.p(1,3) ; 

boot(1,1:3) = [eta beta q] ;%mPart.p(1,1:3) ;   
boot(1,4) = n * log(2 - q) + n * log(beta) - n * beta * log(eta) + ...
    (beta - 1) * sum(log(t)) + 1 / (1 - q) * ...
    sum(log(1 - (1 - q) .* (t ./ eta) .^ beta)) ;%mPart.obj ;

u = rand(n,B) ;
% u = rand(500,B) ;
y = eta * ((1 - (1 - u) .^ ((1 - q) / (2 - q))) ./ (1 - q)) .^ (1 / beta) ;

% KS test (bootstrapped version):
D = zeros(B+1,1) ;
Fn = zeros(n+1,1) ;
for j = 0:n
    Fn(j+1,1) = j / n ;
end
cdf = 1 - (1 - (1 - q) .* (t / eta) .^ beta) .^ ((2 - q) / (1 - q)) ;
cdf = sort(cdf) ;
Dmais = max(abs(Fn(2:n+1) - cdf)) ;
Dmenos = max(abs(Fn(1:n) - cdf)) ;
D(1,1) = max(Dmais, Dmenos) ;

arqBoot = ['./resultado/' arquivo(1:end-4) '/bootAmostras' num2str(B)  '.txt'] ;
dlmwrite(arqBoot, y, 'delimiter', '\t', 'precision', '%.4f') ;
par.nRep = 1 ;
start = tic ;
for i = 1:B 
    % clear part.x part.p part.v part.obj part.evoObj part.mInd part.mIte part.mViz;
    arqResBoot = [arquivo(1:end-4) '/boot' num2str(i) '.txt'] ;
    [mPartBoot] = psoPrincipal(arqResBoot,par,part,y(:,i)) ;
    boot(i+1,:) = [mPartBoot.p(1,1) mPartBoot.p(1,2) ...
        mPartBoot.p(1,3) mPartBoot.obj] ;
    
    % KS test (bootstrapped version):
    cdf = 1 - (1 - (1 - mPartBoot.p(1,3)) .* (y(:,i) / ...
        mPartBoot.p(1,1)) .^ mPartBoot.p(1,2)) .^ ((2 - ...
        mPartBoot.p(1,3)) / (1 - mPartBoot.p(1,3))) ;
    cdf = sort(cdf) ;
    Dmais = max(abs(Fn(2:n+1) - cdf)) ;
    Dmenos = max(abs(Fn(1:n) - cdf)) ;
    D(i+1,1) = max(Dmais, Dmenos) ;
end

elapsedTime = toc(start);
arqBoot2 = ['./resultado/' arquivo(1:end-4) '/bootResultados' num2str(B)  '.txt'] ;
fid = fopen(arqBoot2,'a+') ;
fprintf(fid,'Eta\tBeta\tq\tLog-verossimilhanca\n') ;
fclose(fid)
dlmwrite(arqBoot2, boot, '-append', 'delimiter', '\t', 'precision', '%.4f') ;

minBoot = min(boot) ;
perc5 = quantile(boot,0.05) ;
medianaBoot = median(boot) ;
perc95 = quantile(boot,0.95) ;
maxBoot = max(boot) ;
mediaBoot = mean(boot) ;
desPadBoot = sqrt(var(boot)) ;
compIC95 = perc95 - perc5 ;

arqBoot3 = ['./resultado/' arquivo(1:end-4) '/bootResultados' num2str(B)  '-Analise.txt'] ;
fid = fopen(arqBoot3,'a+') ;
fprintf(fid,'Tempo para rodar %.i amostras bootstrap: %.4f segundos.\n',B,elapsedTime) ;
fprintf(fid,'Eta\tBeta\tq\tLog-verossimilhanca\n') ;
fclose(fid)
dlmwrite(arqBoot3, [minBoot;perc5;medianaBoot;perc95;maxBoot;mediaBoot;desPadBoot], '-append', 'delimiter', '\t', 'precision', '%.4f') ;
fid = fopen(arqBoot3,'a+') ;
fprintf(fid,'\nComprimento dos intervalos\n') ;
fclose(fid) ;
dlmwrite(arqBoot3, compIC95, '-append', 'delimiter', '\t', 'precision', '%.4f') ;

% KS test (bootstrapped version):
perc2p5D = quantile(D,0.025) ;
perc5D = quantile(D,0.05) ;
perc50D = quantile(D,0.50) ;
perc95D = quantile(D,0.95) ;
perc97p5D = quantile(D,0.975) ;
a = find(D > D(1,1)) ;
p = (size(a,1) + 1) / (B + 1) ; 

arqBoot4 = ['./resultado/' arquivo(1:end-4) '/bootResultados' num2str(B)  '-KS.txt'] ;
fid = fopen(arqBoot4,'a+') ;
fprintf(fid,'Teste KS - versao bootstrap:\n') ;
fprintf(fid,'D\tD(p0.025)\tD(p0.05)\tD(p0.5)\tD(p0.95)\tD(p0.975)\tp\n') ;
fclose(fid) ;
dlmwrite(arqBoot4, [D(1,1) perc2p5D perc5D perc50D perc95D perc97p5D p], '-append', 'delimiter', '\t', 'precision', '%.4f') ;

arqBoot5 = ['./resultado/' arquivo(1:end-4) '/bootResultados' num2str(B)  '-KS-distances.txt'] ;
fid = fopen(arqBoot5,'a+') ;
fprintf(fid,'D\n') ;
fclose(fid) ;
dlmwrite(arqBoot5, D, '-append', 'delimiter', '', 'precision', '%.4f') ;
