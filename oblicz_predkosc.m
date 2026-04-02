%% oblicz_predkosc.m
% Skrypt obliczajacy predkosc katowa na podstawie sygnalow enkodera
% inkrementalnego (impulsy_A, impulsy_B).
%
% Metoda: numeryczne przyblizenie pochodnej - wsteczna metoda Eulera
%   omega(k) = (theta(k) - theta(k-1)) / T
%
% Enkoder: 4096 impulsow/obrot, dekodowanie kwadraturowe x4 -> 16384 liczeń/obrot
%
% Czestotliwosci probkowania:
%   f0 = 10 kHz, f1 = 100 kHz, f2 = 1 kHz, f4 = 0.1 kHz

clear; clc; close all;

%% =========================================================
%  PARAMETRY
%% =========================================================
N_PPR    = 4096;           % rozdzielczosc enkodera [impulsy/obrot]
N_COUNTS = 4 * N_PPR;      % dekodowanie x4: liczba zliczen na pelny obrot

% Czestotliwosci probkowania do analizy [Hz]
FS_LIST   = [10e3, 100e3, 1e3, 0.1e3];
FS_LABELS = {'f_0 = 10 kHz', 'f_1 = 100 kHz', 'f_2 = 1 kHz', 'f_4 = 0.1 kHz'};

% Pliki z danymi
PLIKI  = {'prostokat_15rad_s.mat', 'trojkat_20rad_s.mat', 'trojkat_25rad_s.mat'};
TYTULY = {'Prostokatny ~15 rad/s', 'Trojkatny ~20 rad/s', 'Trojkatny ~25 rad/s'};

% Kolory dla kolejnych czestotliwosci
KOLORY = {[0.13 0.47 0.71], [0.84 0.15 0.16], [0.17 0.63 0.17], [0.58 0.40 0.74]};

%% =========================================================
%  TABLICA DEKODOWANIA KWADRATUROWEGO
%  Stan = 2*A + B  (wartosci: 0, 1, 2, 3)
%  Obrot do przodu:  0->2->3->1->0  (+1 za kazde przejscie)
%  Obrot do tylu:    0->1->3->2->0  (-1 za kazde przejscie)
%  quad_lut(stan_poprz+1, stan_curr+1) = zmiana licznika
%% =========================================================
QUAD_LUT = [ 0, -1,  1,  0;   % poprzedni stan = 0 (A=0, B=0)
             1,  0,  0, -1;   % poprzedni stan = 1 (A=0, B=1)
            -1,  0,  0,  1;   % poprzedni stan = 2 (A=1, B=0)
             0,  1, -1,  0];  % poprzedni stan = 3 (A=1, B=1)
LUT_FLAT = QUAD_LUT(:);       % spłaszczona do indeksowania liniowego

%% =========================================================
%  PETLA PO ZBIORACH DANYCH
%% =========================================================
for d = 1:length(PLIKI)
    fprintf('\n===== Przetwarzanie: %s =====\n', PLIKI{d});

    %% ---- Wczytanie danych ----
    dane      = load(PLIKI{d});
    t_raw     = dane.czas(:);
    imp_A     = dane.impulsy_A(:);
    imp_B     = dane.impulsy_B(:);
    omega_odn = dane.omega_odn(:);

    % Czestotliwosc probkowania surowych danych
    fs_raw = 1 / mean(diff(t_raw(1:1000)));
    fprintf('Czestotliwosc surowych danych: %.0f Hz\n', fs_raw);
    fprintf('Zakres omega_odn: [%.2f, %.2f] rad/s\n', ...
            min(omega_odn), max(omega_odn));

    %% ---- Progowanie sygnalow analogowych -> cyfrowe (0/1) ----
    thr_A = 0.5 * (min(imp_A) + max(imp_A));
    thr_B = 0.5 * (min(imp_B) + max(imp_B));
    A_dig = double(imp_A > thr_A);
    B_dig = double(imp_B > thr_B);

    %% ---- Dekodowanie kwadraturowe -> pozycja w zliczeniach ----
    stan    = int32(2*A_dig + B_dig);        % wartosci: 0, 1, 2, 3
    prev_s  = stan(1:end-1) + 1;            % indeks 1..4 (poprzedni)
    curr_s  = stan(2:end)   + 1;            % indeks 1..4 (biezacy)
    idx_lut = (prev_s - 1)*4 + curr_s;      % indeks do LUT_FLAT
    delta   = LUT_FLAT(idx_lut);            % przyrosty zliczen
    pos_raw = [0; cumsum(double(delta))];   % pozycja [zliczenia]

    %% ---- Tworzenie wykresu dla tego zbioru danych ----
    fig = figure('Name', TYTULY{d}, 'NumberTitle', 'off', ...
                 'Position', [30, 30, 1400, 800]);
    sgtitle(sprintf('Predkosc katowa z enkodera  —  Sygnal %s', TYTULY{d}), ...
            'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

    %% ---- Petla po czestotliwosciach probkowania ----
    for f = 1:length(FS_LIST)
        fs_t = FS_LIST(f);
        T    = 1 / fs_t;          % okres probkowania [s]

        % Wektor czasu dla docelowej czestotliwosci
        t_target = (0 : floor(t_raw(end) * fs_t))' / fs_t;
        % Przycinamy do zakresu dostepnych danych
        t_target = t_target(t_target <= t_raw(end));

        % --- Probkowanie pozycji metoda ZOH (zero-order hold) ---
        % Symuluje odczyt licznika enkodera w chwilach probkowania
        pos_t = interp1(t_raw, pos_raw, t_target, 'previous', 'extrap');

        % --- Wsteczna metoda Eulera ---
        % omega(k) = [theta(k) - theta(k-1)] / T
        % theta [rad] = pos [zliczenia] * (2*pi / N_COUNTS)
        d_pos = [0; diff(pos_t)];
        omega_enc = d_pos * (2 * pi / N_COUNTS) / T;   % [rad/s]

        % --- Interpolacja przebiegu referencyjnego do tej samej siatki ---
        omega_ref_t = interp1(t_raw, omega_odn, t_target, 'linear', 'extrap');

        % --- Subplot ---
        subplot(2, 2, f);
        plot(t_target, omega_ref_t, 'k--', 'LineWidth', 1.8, ...
             'DisplayName', '\omega_{odn}');
        hold on;
        plot(t_target, omega_enc, 'Color', KOLORY{f}, 'LineWidth', 1.0, ...
             'DisplayName', ['\omega_{enc},\ ', FS_LABELS{f}]);
        hold off;

        xlabel('Czas [s]',     'FontSize', 10);
        ylabel('\omega [rad/s]', 'FontSize', 10);
        title(FS_LABELS{f}, 'FontSize', 11, 'FontWeight', 'bold', ...
              'Interpreter', 'tex');
        legend('show', 'Location', 'best', 'FontSize', 9, 'Interpreter', 'tex');
        grid on;
        xlim([0, t_target(end)]);

        fprintf('  fs = %8g Hz  |  N_probek = %d\n', fs_t, length(t_target));
    end

    % Zapis wykresu jako PNG
    nazwa_out = sprintf('wykresy_%s', strrep(PLIKI{d}, '.mat', ''));
    saveas(fig, [nazwa_out '.png']);
    fprintf('Zapisano: %s.png\n', nazwa_out);
end

fprintf('\n--- Obliczenia zakonczone. ---\n');
