% Parâmetros canal real
medicao = load('Prx_Real_2021_1.mat');
% Extraindo os dados e transpondo os vetores
vtPrxdBm = medicao.dPrx'; 
vtDist = medicao.dPath';
% Tamanho da janela
dW = 5;
% Transforma potência em mWatts
vtPtrxmW = 10.^(vtPrxdBm/10);
nSamples = length(vtPtrxmW);
% Vetores para canal estimado
vtDesLarga = [];
vtDesPequeEst = [];
%
% Cálculo do desvanecimenro lento e rápido
dMeiaJanela = round((dW-1)/2);  % Meia janela
ij = 1;
for ik = dMeiaJanela + 1 : nSamples - dMeiaJanela
    % Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
    vtDesLarga(ij) = 10*log10(mean(vtPtrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
    % Desvanecimento de pequena escala [dB]
    vtDesPequeEst(ij) = vtPrxdBm(ik)-vtDesLarga(ij);
    ij = ij + 1;
end
%
% Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
indexes = dMeiaJanela+1 : nSamples-dMeiaJanela;
%vtPrxW = ((10.^(vtPrxdBm(indexes)./10))/1000);
vtPtrxmWNew = 10.^(vtPrxdBm(indexes)/10);
desLarga_Lin = (10.^(vtDesLarga(1:length(indexes))./10));
envNormal = sqrt(vtPtrxmWNew)./sqrt(desLarga_Lin);
%
% Ajuste no tamanho dos vetores devido a filtragem
vtDistEst = vtDist( dMeiaJanela+1 : nSamples-dMeiaJanela );
vtPrxdBm = vtPrxdBm( dMeiaJanela+1 : nSamples-dMeiaJanela );
%
% Cálculo reta de perda de percurso
vtDistLog = log10(vtDist);
vtDistLogEst = log10(vtDistEst);
% Cálculo do coeficientes da reta que melhor se caracteriza a perda de percurso
dCoefReta = polyfit(vtDistLogEst,vtPrxdBm,1); 
% Expoente de perda de percurso estimado
dNEst = -dCoefReta(1)/10;
disp(['Estimação dos parâmetros de larga escala (W = ' num2str(dW) '):'])
disp(['   Expoente de perda de percurso estimado n = ' num2str(dNEst)]);
% Perda de percurso estimada para os pontos de medição
vtPathLossEst = polyval(dCoefReta,vtDistLogEst);  
%
% Sombreamento
vtShadCorrEst = vtDesLarga - vtPathLossEst;
% Calcula a variância do sombreamento estimado
stdShad = std(vtShadCorrEst);
meanShad = mean(vtShadCorrEst);
disp(['   Desvio padrão do sombreamento estimado = ' num2str(stdShad)]);
disp(['   Média do sombreamento estimado = ' num2str(meanShad)]);
vtPathLossEst = - vtPathLossEst;
%vtPrxEst = sPar.txPower - vtPathLossEst + vtShadCorrEst + vtDesPequeEst;
% vtDesLarga - vtPathLossEst - vtPathLossEst
% Figuras do canal estimado
figure;
% Potência recebida com canal completo
plot(vtDistLogEst,vtPrxdBm); hold all;
% Potência recebida com path loss
plot(vtDistLogEst,-vtPathLossEst,'linewidth', 2)
% Potência recebida com path loss e shadowing
plot(vtDistLogEst,-vtPathLossEst+vtShadCorrEst,'linewidth', 2)
%title('Canal estimado: Potência recebida no receptor vs. log da distância')
xlabel('log_{10}(d)');
ylabel('Potência [dBm]');
legend('Prx canal completo', 'Prx (somente perda de percurso)', 'Prx (perda de percurso + sombreamento)');
title('Prx Medição Real vs Estimada');
%

