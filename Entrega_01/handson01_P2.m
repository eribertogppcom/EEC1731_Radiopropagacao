close all;clear all;clc;
addpath('./CODES/HD_03_MATLAB/fitmethis')
% Parâmetros canal real
medicao = load('Prx_Real_2021_1.mat');
% Extraindo os dados e transpondo os vetores
vtPrxdBm = medicao.dPrx'; 
vtDist = medicao.dPath';
%
% Várias janelas de filtragem para testar a estimação
vtW = [2 5 10];
dNEst = [];
stdShad =[];
meanShad = [];

for iw = 1: length(vtW)
    % Configura valor da janela de filtragem
    dW = vtW(iw);
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
dNEst(iw) = -dCoefReta(1)/10;
% Perda de percurso estimada para os pontos de medição
vtPathLossEst = polyval(dCoefReta,vtDistLogEst);  
% Sombreamento
vtShadCorrEst = vtDesLarga - vtPathLossEst;
% Calcula a variância do sombreamento estimado
stdShad(iw) = std(vtShadCorrEst);
meanShad(iw) = mean(vtShadCorrEst);
vtPathLossEst = - vtPathLossEst;
end

% Estimação cega via FitMeThis
disp(' ')
disp('Estimação do Fading para várias janelas (estudo númerico sem conhecimento a priori do canal)');
disp('Resultados com FITMETHIS')
for iw = 1:length(vtW)%
    disp(['Janela W = ' num2str(vtW(iw))]);%Tabulando o resultado
    %fprintf('    Janela     Primeira melhor PDF   Parâmetro(s) da primeira melhor PDF 	Segunda melhor PDF 	Parâmetro(s) da segunda melhor PDF\n');
    fprintf('------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf('Estimação dos parâmetros de larga escala (W):%3d\nExpoente de perda de percurso estimado (n):%3d\nDesvio padrão do sombreamento estimado :%3d\nMédia do sombreamento estimado :%3d\n', vtW(iw), dNEst(iw), stdShad(iw), meanShad(iw));       
    %
    fprintf('------------------------------------------------------------------------------------------------------------------------------\n');
    fitWindow{iw} = fitmethis([envNormal]','figure','off');
    
    disp(' ')
end
