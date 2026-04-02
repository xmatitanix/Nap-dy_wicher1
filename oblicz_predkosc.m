%% oblicz_predkosc.m
% =========================================================================
%  CYFROWY POMIAR PREDKOSCI KATOWEJ Z ENKODERA INKREMENTALNEGO
% =========================================================================
%  Implementuje dwie metody i porownuje je z sygnalem referencyjnym:
%    (M)   Metoda M  — Euler wsteczny (staly okres, zliczanie impulsow)
%    (M/T) Metoda M/T zmodyfikowana  — staly okres, interpolacja podprobkowa
%
%  Literatura:
%   [1] S. Brock, K. Zawirski: "Cyfrowy pomiar predkosci obrotowej
%       w napedzie elektrycznym", PAR 1/2005
%   [2] S. Brock, J. Deskur: "The problem of measurement and control
%       of speed in a drive with an inaccurate measuring position
%       transducer", IECON 2008
%
% =========================================================================
%  TEORIA
% =========================================================================
%
%  1. ENKODER INKREMENTALNY — DEKODOWANIE KWADRATUROWE 4x
%  -------------------------------------------------------
%  Enkoder generuje dwa sygnaly kwadratury A i B, przesuniete w fazie o 90°.
%  Dekodowanie 4x (oba zbocza obu kanalow) daje rozdzielczosc:
%
%    N_COUNTS = 4 * N_PPR   [zliczen/obrot]                         (1)
%
%  Tablica stanow: stan = 2*A + B in {0,1,2,3}
%    Obrot CW (do przodu): 0->2->3->1->0  => +1 na kazde przejscie
%    Obrot CCW (do tylu):  0->1->3->2->0  => -1 na kazde przejscie
%
%  Polozenie katowe w radianach:
%    theta(k) = count(k) * (2*pi / N_COUNTS)   [rad]                (2)
%
%  2. METODA M — EULER WSTECZNY (staly okres probkowania Tp)
%  ----------------------------------------------------------
%  Numeryczne przyblizenie pochodnej polozenia po czasie:
%
%    omega_M(k) = [theta(k) - theta(k-1)] / Tp                      (3)
%               = Delta_count(k) * (2*pi / N_COUNTS) / Tp  [rad/s]
%
%  Blad kwantyzacji: Delta_count(k) to liczba calkowita => blad +/-1.
%  Wzgledny blad rozdzielczosci:
%
%    epsilon_M = 1 / |Delta_count_sr|                               (4)
%
%  gdzie Delta_count_sr = omega * N_COUNTS * Tp / (2*pi).
%  Wniosek: przy duzej predkosci lub dlugim Tp blad jest maly.
%  Przy malej predkosci (Delta_count -> 0) blad staje sie dominujacy.
%
%  3. METODA M/T ZMODYFIKOWANA (Brock & Zawirski 2005, rown. 4–5 w [1])
%  ----------------------------------------------------------------------
%  Idea: korekcja granicy okna czasowego z dokladnoscia do czasu miedzy
%  kolejnymi impulsami enkodera (interpolacja podprobkowa).
%
%  Licznik L2 (w oprogramowaniu: zmienna C) liczy takty zegara f_g od
%  ostatniego impulsu enkodera. Przy kazdym impulsie jest zerowany.
%
%  Rzeczywisty czas pomiaru:
%    Tp_actual(k) = Tp + [C(k-1) - C(k)]                            (5)
%
%  gdzie C(k) = czas od ostatniego impulsu enkodera do k-tej chwili
%  probkowania [w sekundach lub w taktach zegara].
%
%  Predkosc:
%    omega_MT(k) = Delta_count(k) * (2*pi / N_COUNTS) / Tp_actual(k) (6)
%
%  Metoda zapewnia:
%    - staly nominalny okres probkowania (jak M)
%    - lepsza dokladnosc statyczna dzieki interpolacji czasu (jak T)
%    - poprawne dzialanie przy predkosci rownej zeru

clear; clc; close all;

% =========================================================================
%  PARAMETRY
% =========================================================================
N_PPR    = 4096;           % rozdzielczosc enkodera [impulsy/obrot]
N_COUNTS = 4 * N_PPR;      % zliczenia x4 kwadratura = 16384 [zliczen/obrot]

% Czestotliwosci probkowania do analizy — wg polecenia: f0..f4
FS_LIST   = [10e3,    100e3,    1e3,    0.1e3];
FS_LABELS = {'f_0 = 10 kHz', 'f_1 = 100 kHz', 'f_2 = 1 kHz', 'f_4 = 0.1 kHz'};

% Pliki wejsciowe
PLIKI  = {'prostokat_15rad_s.mat', 'trojkat_20rad_s.mat', 'trojkat_25rad_s.mat'};
TYTULY = {'Prostokatny ~15 rad/s', 'Trojkatny ~20 rad/s', 'Trojkatny ~25 rad/s'};

% Paleta kolorow dla czterech czestotliwosci
KOLORY_M  = {[0.13 0.47 0.71], [0.84 0.15 0.16], ...
             [0.17 0.63 0.17], [0.58 0.40 0.74]};
KOLORY_MT = {[0.05 0.25 0.50], [0.60 0.05 0.05], ...
             [0.05 0.40 0.05], [0.35 0.15 0.55]};

% =========================================================================
%  TABLICA DEKODOWANIA KWADRATUROWEGO
%  QUAD_LUT(wiersz = stan_poprzedni+1, kolumna = stan_biezacy+1)
%  Uzycie sub2ind gwarantuje poprawne row-major mapowanie w MATLAB
%  (MATLAB(:) splaszcza kolumnowo, wiec reczna arytmetyka daje odwrocone znaki!)
% =========================================================================
QUAD_LUT = [ 0, -1,  1,  0;   % stan_poprz = 0: A=0,B=0
             1,  0,  0, -1;   % stan_poprz = 1: A=0,B=1
            -1,  0,  0,  1;   % stan_poprz = 2: A=1,B=0
             0,  1, -1,  0];  % stan_poprz = 3: A=1,B=1

% =========================================================================
%  PETLA PO ZBIORACH DANYCH
% =========================================================================
for d = 1:length(PLIKI)
    fprintf('\n===== Plik: %s =====\n', PLIKI{d});
    tic;

    % ---------------------------------------------------------------------
    %  Wczytanie danych
    % ---------------------------------------------------------------------
    dane      = load(PLIKI{d});
    t_raw     = dane.czas(:);
    imp_A     = dane.impulsy_A(:);
    imp_B     = dane.impulsy_B(:);
    omega_odn = dane.omega_odn(:);

    % Czestotliwosc probkowania surowych danych (oczekiwana ~250 kHz)
    fs_raw = 1.0 / mean(diff(t_raw(1:1000)));
    fprintf('  fs_raw = %.0f Hz  |  T_raw = %.2f us\n', ...
            fs_raw, 1e6/fs_raw);
    fprintf('  Zakres omega_odn: [%.2f, %.2f] rad/s\n', ...
            min(omega_odn), max(omega_odn));

    % ---------------------------------------------------------------------
    %  Progowanie sygnalu analogowego -> cyfrowy 0/1
    %  Prog = polowa zakresu (metoda polsumowania min i max)
    %  Odporne na wolnozmienny dryft DC sygnalu
    % ---------------------------------------------------------------------
    thr_A = 0.5 * (min(imp_A) + max(imp_A));
    thr_B = 0.5 * (min(imp_B) + max(imp_B));
    A_dig = double(imp_A > thr_A);   % sygnal cyfrowy kanalA: 0 lub 1
    B_dig = double(imp_B > thr_B);   % sygnal cyfrowy kanalB: 0 lub 1

    % ---------------------------------------------------------------------
    %  Dekodowanie kwadraturowe 4x — w pelni wektorowe (brak petli)
    %  Zlozonosc: O(N), gdzie N = liczba probek surowych (~10^6)
    % ---------------------------------------------------------------------
    stan   = 2*A_dig + B_dig;                    % stan in {0,1,2,3}
    prev_s = stan(1:end-1) + 1;                  % poprzedni stan (1-indeks)
    curr_s = stan(2:end)   + 1;                  % biezacy stan   (1-indeks)
    % sub2ind: poprawny dostep QUAD_LUT(poprz, biezacy) — patrz komentarz LUT
    delta  = QUAD_LUT(sub2ind(size(QUAD_LUT), prev_s, curr_s));
    % Pozycja = skumulowana suma zmian zliczen od probki 0
    pos_raw = [0; cumsum(delta)];                % [zliczenia], dlugosc N

    % ---------------------------------------------------------------------
    %  Preobliczenie czasow przejsc (do metody M/T)
    %  Przejscie = chwila, w ktorej licznik kwadraturowy sie zmienil
    %  Uzywane do wyznaczenia C(k) — czasu od ostatniego impulsu enkodera
    % ---------------------------------------------------------------------
    maska_przejsc = [false; delta ~= 0];         % logiczna maska przejsc
    t_trans       = t_raw(maska_przejsc);        % czasy przejsc [s]
    % Zapewnienie monotonii (usuniecie ewentualnych duplikatow)
    t_trans       = unique(t_trans);

    fprintf('  Liczba przejsc kwadraturowych: %d (sr. %.0f/s)\n', ...
            length(t_trans), length(t_trans)/t_raw(end));
    fprintf('  Czas przetwarzania wstepnego: %.2f s\n', toc);

    % =====================================================================
    %  TWORZENIE FIGURY
    % =====================================================================
    fig = figure('Name', TYTULY{d}, 'NumberTitle', 'off', ...
                 'Position', [30, 30, 1400, 900]);
    sgtitle(sprintf('Predkosc katowa z enkodera — %s', TYTULY{d}), ...
            'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

    % =====================================================================
    %  PETLA PO CZESTOTLIWOSCIACH PROBKOWANIA
    % =====================================================================
    for f = 1:length(FS_LIST)
        fs_t = FS_LIST(f);
        Tp   = 1.0 / fs_t;   % nominalny okres probkowania [s]

        % -----------------------------------------------------------------
        %  Wektor czasu docelowego
        % -----------------------------------------------------------------
        t_target = (0 : floor(t_raw(end) * fs_t))' / fs_t;
        t_target = t_target(t_target <= t_raw(end));
        N_t      = length(t_target);

        % -----------------------------------------------------------------
        %  ZOH: odczyt polozenia licznika w chwilach probkowania
        %  floor(t * fs_raw) + 1 to indeks ostatniej probki surowej <= t
        %  Odpowiednik odczytu sprzetowego licznika impulsow
        % -----------------------------------------------------------------
        idx_zoh = min(floor(t_target * fs_raw) + 1, length(pos_raw));
        pos_t   = pos_raw(idx_zoh);           % pozycja [zliczenia] w t_target
        d_pos   = [0; diff(pos_t)];           % Delta_count(k) wg rownania (3)

        % -----------------------------------------------------------------
        %  METODA M — EULER WSTECZNY  [rown. 3]
        %  omega_M(k) = Delta_count(k) * (2*pi/N_COUNTS) / Tp
        % -----------------------------------------------------------------
        omega_M = d_pos * (2*pi / N_COUNTS) / Tp;   % [rad/s]

        % -----------------------------------------------------------------
        %  METODA M/T ZMODYFIKOWANA  [rown. 5–6]
        %
        %  Krok 1: wyznacz C(k) — czas od ostatniego impulsu enkodera
        %          do k-tej chwili probkowania [s]
        %  interp1 z metoda 'previous': dla kazdego t_target(k) zwraca
        %  ostatni t_trans(i) <= t_target(k)  — operacja O(N_t * log(M))
        % -----------------------------------------------------------------
        t_last = interp1(t_trans, t_trans, t_target, 'previous');
        % Dla probek przed pierwszym przejsciem: t_last = NaN -> C=0
        t_last(isnan(t_last)) = t_target(isnan(t_last));
        C_k    = t_target - t_last;           % C(k) >= 0  [s]

        %  Krok 2: Tp_actual(k) = Tp + C(k-1) - C(k)  [rown. 5]
        C_prev     = [C_k(1); C_k(1:end-1)];
        Tp_actual  = Tp + C_prev - C_k;      % rzeczywisty czas pomiaru [s]

        %  Krok 3: omega_MT(k) = Delta_count(k) * (2*pi/N_COUNTS) / Tp_actual(k)
        %  Zabezpieczenie: jesli Tp_actual < Tp/4, podstaw 0
        %  (sytuacja gdy blad podprobkowy przekracza 75% okresu probkowania)
        omega_MT              = d_pos * (2*pi / N_COUNTS) ./ Tp_actual;
        omega_MT(Tp_actual < Tp/4) = 0;      % [rad/s]

        % -----------------------------------------------------------------
        %  Interpolacja sygnalu referencyjnego omega_odn do siatki docelowej
        % -----------------------------------------------------------------
        omega_ref_t = interp1(t_raw, omega_odn, t_target, 'linear', 'extrap');

        % -----------------------------------------------------------------
        %  Wyswietlanie wynikow
        % -----------------------------------------------------------------
        subplot(2, 2, f);

        plot(t_target, omega_ref_t, 'k--', 'LineWidth', 1.8, ...
             'DisplayName', '\omega_{odn}');
        hold on;
        plot(t_target, omega_M,  'Color', KOLORY_M{f},  'LineWidth', 0.9, ...
             'DisplayName', ['Metoda M: ', FS_LABELS{f}]);
        plot(t_target, omega_MT, 'Color', KOLORY_MT{f}, 'LineWidth', 1.3, ...
             'LineStyle', '-.', ...
             'DisplayName', ['Metoda M/T: ', FS_LABELS{f}]);
        hold off;

        xlabel('Czas [s]',       'FontSize', 10);
        ylabel('\omega [rad/s]', 'FontSize', 10);
        title(FS_LABELS{f}, 'FontSize', 11, 'FontWeight', 'bold', ...
              'Interpreter', 'tex');
        legend('Location', 'best', 'FontSize', 8, 'Interpreter', 'tex');
        grid on;
        xlim([0, t_target(end)]);

        fprintf('  [%s] N=%5d | M: [%6.1f, %6.1f] | M/T: [%6.1f, %6.1f] rad/s\n', ...
                FS_LABELS{f}, N_t, min(omega_M), max(omega_M), ...
                min(omega_MT), max(omega_MT));
    end

    % Zapis wykresu
    out_png = sprintf('wykresy_%s.png', strrep(PLIKI{d}, '.mat', ''));
    saveas(fig, out_png);
    fprintf('  Zapisano: %s  (czas calkowity: %.2f s)\n', out_png, toc);
end

fprintf('\n--- Obliczenia zakonczone. ---\n');
