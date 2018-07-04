%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                    Experiencia 4                    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definindo constantes:

f_a = 100e6; % Taxa de amostragem do sinal
T_a = 1/(f_a); % Periodo de amostragem correspondente
f_0 = 10e6; % Frequencia central do pulso senoidal transmitido

% ... do item 1
TAU = 10e-6; % Largura do pulso senoidal transmitido

% ... do item 3
T = 200e-6; % Largura da janela observada
T_d = 20e-6; % Tempo de atraso do sinal recebido
A = 1.64365; % Amplitude do sinal recebido
sigma = 0.1; % Desvio-padraao do ruiido gerado

% ... do item 5
BETA = 0.1; % Paraametro da funccaao usada para gerar pulso chirp; valor arbitraario

% constantes auxiliares
c = 3e8; % Velocidade da luz, [m/s]
R = (c * T_d)/2; % Distancia do alvo ao radar, no momento em que o pulso incide no alvo

vet_t_am = 0:T_a:TAU; % Vetor com os tempos nos quais ee realizada uma amostragem

%%%%%%%%%% Parte Experimental 1 - Usando pulso senoidal %%%%%%%%%%%%%%%
%% 1 - Gerando pulso senoidal:
p_1 = sin(2*pi*f_0*vet_t_am);

%% 2 - Gerando filtro casado:
h_n = conj(p_1(end:-1:1));

%% 3 - Gerando o sinal recebido
% Adequando o sinal sem ruiido
vet_atraso = zeros(1, int32(T_d/T_a)); % amostras nulas antes da chegada do sinal
p_1_atrasado = [vet_atraso, p_1]; % sinal com o atraso, de tamanho naao fixo
tamanho_janela = int32(T/T_a);
resto_da_janela = zeros(1, (tamanho_janela - size(p_1_atrasado, 2)));
p_1_atr_janelado = [p_1_atrasado, resto_da_janela];

% Gerando e adicionando ruiido gaussiano
ruido_gaussiano = sigma.*randn(1, tamanho_janela);
ruido_gaussiano_curto = sigma.*randn(1, size(p_1, 2));

% Caalculo da poteencia do sinal (com amplitude da senooide igual a 1)
pot_sinal_puro = rms(p_1)^2;
pot_ruido = rms(ruido_gaussiano)^2;

% Caalculo da amplitude A tal que as poteencias meedias do sinal e do ruiido
% sejam iguais - situaccaao que equivale a SNR_dB = 0[dB]
A = sqrt(pot_ruido/pot_sinal_puro); % Observe que a poteencia do ruiido
                                    % seraa igual aa variaancia, em meedia

% Gerando sinal com amplitude corrigida
p_1_atr_jan_amp = A*p_1_atr_janelado;

% Verificando que, de fato, este valor de A cria um sinal que obedece o
% enunciado
relacao_s_r = snr((A*p_1), ruido_gaussiano_curto);
fprintf('Item 3 (senóide): Verificando a relação sinal-ruído\n');
fprintf('SNR [dB]: %.2f\n\n', relacao_s_r);

% Por fim, gerando o sinal recebido, x_n
x_n = (A*p_1_atr_janelado) + ruido_gaussiano;

%% 4 - Comparaccooes
% Calculando saiidas do filtro casado
saida_amp_corrigida = filter(h_n, 1, p_1_atr_jan_amp);
saida_x_n = filter(h_n, 1, x_n);

% Visualizaccooes
figure(1);
plot(saida_amp_corrigida, 'r')
title('Saída do filtro (senóide; entrada sem ruído; amplitude corrigida)');

figure(2);
plot(saida_x_n, 'g');
title('Saída do filtro (senóide; entrada ruidosa; amplitude corrigida)');

x = 1:1:size(saida_x_n, 2);

figure(3);
plot(x, saida_x_n, 'g', x, saida_amp_corrigida, 'r');
title('Comparação de resultados (senóide)');

% Calculando as posiccooes e tempos do maaximo e do atraso
n_0 = double((find(ismember(saida_x_n, max(saida_x_n(:)))) - 1));
t_0 = (n_0/size(x, 2))*T;

n_atraso = T_d/T_a;
t_atraso = (n_atraso/size(x, 2))*double(T);

delta_n = abs(n_atraso - n_0);
delta_t = abs(t_atraso - t_0);
fprintf('Item 4 (senóide): Verificando diferenças nos tempos de detecção\n');
fprintf('delta_n: %d  [amostras]\n', int32(delta_n));
fprintf('delta_t: %.2f [microssegundos]\n\n', delta_t*10^6);

%% 5 - Repetindo passos acima para chirp
%%%%%%%%%%%%%%% Parte Experimental 2 - Usando chirp %%%%%%%%%%%%%%%%%
%% 1' - Gerando chirp:
p_2 = zeros(1, size(p_1, 2));
for n = 1:size(p_2, 2)
    p_2(n) = sin(2*pi*f_0*(1 + (BETA/TAU)*n)*n);
end

%% 2' - Gerando filtro casado:
h_2_n = conj(p_2(end:-1:1));

%% 3' - Gerando o sinal recebido
% Adequando o sinal sem ruiido
p_2_atrasado = [vet_atraso, p_2]; % sinal com o atraso, de tamanho naao fixo
tamanho_janela = int32(T/T_a);
resto_da_janela = zeros(1, (tamanho_janela - size(p_2_atrasado, 2)));
p_2_atr_janelado = [p_2_atrasado, resto_da_janela];

% Gerando e adicionando ruiido gaussiano
ruido_gaussiano = sigma.*randn(1, tamanho_janela);
ruido_gaussiano_curto = sigma.*randn(1, size(p_2, 2));

% Caalculo da poteencia do sinal (com amplitude da senooide igual a 1)
pot_sinal_puro_2 = rms(p_2)^2;

% Caalculo da amplitude A tal que as poteencias meedias do sinal e do ruiido
% sejam iguais - situaccaao que equivale a SNR_dB = 0[dB]
A = sqrt(pot_ruido/pot_sinal_puro_2); % Observe que a poteencia do ruiido
                                    % seraa igual aa variaancia, em meedia

% Gerando sinal com amplitude corrigida
p_2_atr_jan_amp = A*p_2_atr_janelado;

% Verificando que, de fato, este valor de A cria um sinal que obedece o
% enunciado
relacao_s_r_2 = snr((A*p_2), ruido_gaussiano_curto);
fprintf('Item 3 (chirp): Verificando a relação sinal-ruído\n');
fprintf('SNR [dB]: %.2f\n\n', relacao_s_r_2);

% Por fim, gerando o sinal recebido, x_2_n
x_2_n = (A*p_2_atr_janelado) + ruido_gaussiano;

%% 4' - Comparaccooes
% calculando saiidas do filtro casado
saida_amp_corrigida_2 = filter(h_2_n, 1, p_2_atr_jan_amp);
saida_x_2_n = filter(h_2_n, 1, x_2_n);

% Visualizaccooes
figure(4);
plot(saida_amp_corrigida_2, 'r')
title('Saída do filtro (chirp; entrada sem ruído; amplitude corrigida)');

figure(5);
plot(saida_x_2_n, 'g');
title('Saída do filtro (chirp; entrada ruidosa; amplitude corrigida)');

x = 1:1:size(saida_x_2_n, 2);

figure(6);
plot(x, saida_x_2_n, 'g', x, saida_amp_corrigida_2, 'r');
title('Comparação de resultados (chirp)');

% Calculando as posiccooes e tempos do maaximo e do atraso
n2_0 = double((find(ismember(saida_x_2_n, max(saida_x_2_n(:)))) - 1));
t2_0 = (n2_0/size(x, 2))*T;

t2_atraso = (n_atraso/size(x, 2))*double(T);

delta_n2 = abs(n_atraso - n2_0);
delta_t2 = abs(t2_atraso - t2_0);
fprintf('Item 4 (chirp): Verificando diferenças nos tempos de detecção\n');
fprintf('delta_n2: %d  [amostras]\n', int32(delta_n2));
fprintf('delta_t2: %.2f [microssegundos]\n\n', delta_t2*10^6);
