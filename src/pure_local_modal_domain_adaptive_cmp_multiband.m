

%% ============================================================
% local_modal_multiband_pure_local_modal_adaptive_ceil_uniform_cmp_phasecheck.m
%
% [순수 로컬 모달 + kR ceil 기반 주파수별 적응 모달 차수 + 균일 CMP 출력 제약 버전]
%
% 목적:
%   1) 스피커 위치 선택은 시뮬레이션 모델 기반으로 먼저 수행
%   2) 목적함수는 끝까지 local modal coefficient matching 유지
%   3) 이전 측정에서 4채널 상대 응답이 동일하게 확인되었으므로,
%      외부 qSafe_freq 파일을 불러오지 않고 균일 CMP 출력 제약을 직접 사용
%   4) 이후 실제 설치 후 boundary H를 측정해서 이 결과와 비교 가능
%   5) 고정 모달 차수 대신 각 주파수의 kR 상계(ceil(kR))와 Bessel 안정성에 따라
%      사용할 modal index를 자동 선택한다.
%
% 핵심 구조:
%   - target pressure on boundary -> local modal coefficient beta
%   - candidate pressure response on boundary -> local modal response B
%   - selection: multi-frequency modal-domain greedy LS
%   - final q(f): modal-domain LS + uniform CMP power constraint
%
% 필요한 외부 qSafe 파일 없음.
%   - 이전 측정에서 qSafe_freq(:,f)가 모든 채널에서 동일하게 나왔으므로,
%     qSafe_uniform_value = 10.5를 직접 사용한다.
%   - 이는 기존 qSafe_final_scale=30, qSafe_freq=0.35일 때의
%     effective qSafe = 10.5와 동일하다.
%   - apply_weighted_constraints 기준으로 ||q||^2 <= 10.5^2 = 110.25
%     수준의 균일 CMP-style 출력 제한과 같다.
%   - 기존 unwrap phase plot은 유지하고, wrapped phase error plot을 추가로 출력한다.
%   - 1500 Hz에서 gain-scaled boundary/SPL diagnostic을 추가로 출력한다.
%% ============================================================

%% [1] 기본 물리 세팅
c = 343;

room.Lx = 9.2;
room.Ly = 6.5;

% 방 평균 흡음계수 가정
room.alpha_wall = 0.7;

% 압력 반사계수. alpha_wall이 흡음률이면 반사 진폭은 sqrt(1-alpha)
room.beta_wall = sqrt(1 - room.alpha_wall);

% [왼쪽 벽, 오른쪽 벽, 아래쪽 벽, 위쪽 벽]
room.refl = room.beta_wall * [1, 1, 1, 1];

%% [2] 구역 설정
bright.center = [4.6, 2.75];
bright.radius_outer = 0.10;
bright.radius_inner = 0.05;

dark.center = [4.6, 5.9];
dark.radius = 0.10;

%% [3] 스피커 후보군 설정
array.center = bright.center;
array.radius = 1.2;
anglesDeg = 30:-4:-210;

nSelect = 4;

%% [4] 알고리즘 파라미터
lambda = 1e-5;
darkPenalty = 2.5;

% local modal block 가중치
w_outer = 1.0;
w_inner = 0.6;

% 주파수 설정
fList = 100:10:2000;
K = numel(fList);
freqWeight = ones(K,1);

% selection 단계에서는 아직 선택된 4채널이 확정되기 전이므로
% 매우 느슨한 균일 제약을 사용한다.
qSafe_selection_default = 100.0;

% final q 단계에서는 외부 relative_acoustic_constraint_result.mat를 불러오지 않는다.
% 이전 bright-center IR 측정에서 4채널 qSafe가 모두 동일하게 확인되었으므로,
% 채널별 차등 제약이 아니라 균일 CMP-style 출력 제한으로 정리한다.
% 기존 설정: qSafe_final_scale = 30, qSafe_freq = 0.35
% 따라서 effective qSafe = 30 * 0.35 = 10.5
qSafe_uniform_value = 10.5;
qSafe_equivalent_Pmax = qSafe_uniform_value^2;

% 모든 선택 채널과 모든 주파수에 동일한 effective qSafe 적용
qSafe_effective = qSafe_uniform_value * ones(nSelect, K);

% 시각화용 격자
Nx = 220;
Ny = 220;

%% [5] 후보 스피커 생성
% 모든 후보 스피커는 bright center를 향함
candidates = generate_arc_candidates(array.center, array.radius, anglesDeg);
Nc = numel(candidates);

%% [6] boundary 샘플링
N_boundary = 16;
th = linspace(0, 2*pi, N_boundary + 1).';
th(end) = [];

brightOuterPts = [bright.center(1) + bright.radius_outer*cos(th), ...
                  bright.center(2) + bright.radius_outer*sin(th)];

brightInnerPts = [bright.center(1) + bright.radius_inner*cos(th), ...
                  bright.center(2) + bright.radius_inner*sin(th)];

darkBoundaryPts = [dark.center(1) + dark.radius*cos(th), ...
                   dark.center(2) + dark.radius*sin(th)];

brightCenterPt = bright.center;
darkCenterPt = dark.center;

%% [7] 로컬 모달 차수 설정: 주파수별 adaptive N
% 기존 고정값: outer/dark M=4, inner M=2.
% 문제점: 낮은 주파수 또는 Bessel 함수가 작은 주파수에서 1/J_m(kR) 항이 과도하게 커질 수 있음.
% 수정: 각 주파수마다 M(f,R)=ceil(kR)을 기본 상계로 두고, |J_m(kR)|를 보고 사용할 mode index를 자동 선택.
%
% Mmax_*는 허용 가능한 최대 차수이고, 실제 사용 차수/인덱스는 주파수별로 달라진다.
Mmax_outer = 4;
Mmax_inner = 2;
Mmax_dark  = 4;

modalOrder.besselAbsMin = 2e-2;   % |J_m(kR)|가 너무 작은 mode는 제외
modalOrder.regEta = 1e-6;         % modal transform LS regularization

numModes_outer = zeros(K,1);
numModes_inner = zeros(K,1);
numModes_dark  = zeros(K,1);

%% [8] 목표 음장 설정
% bright center 기준 x- / x+ 방향 0.6 m 위치에 이상적 점음원 2개
vs_offset = 0.60;
virtual_source_L = [bright.center(1) - vs_offset, bright.center(2)];
virtual_source_R = [bright.center(1) + vs_offset, bright.center(2)];

%% [9] 주파수별 데이터 저장 변수
beta_total_all = cell(K,1);
B_total_all = cell(K,1);
targetData = cell(K,1);
modalInfo = cell(K,1);

fprintf('=== 주파수별 local modal target 및 candidate block 생성 시작 ===\n');

%% [10] 주파수별 target beta와 candidate modal response 생성
for kk = 1:K
    f = fList(kk);
    k = 2*pi*f/c;

    m_idx_outer = choose_modal_indices(bright.radius_outer, k, Mmax_outer, modalOrder.besselAbsMin);
    m_idx_inner = choose_modal_indices(bright.radius_inner, k, Mmax_inner, modalOrder.besselAbsMin);
    m_idx_dark  = choose_modal_indices(dark.radius,        k, Mmax_dark,  modalOrder.besselAbsMin);

    numModes_outer(kk) = length(m_idx_outer);
    numModes_inner(kk) = length(m_idx_inner);
    numModes_dark(kk)  = length(m_idx_dark);

    F_outer = local_modal_transform_regls(bright.radius_outer, th, k, m_idx_outer, modalOrder.regEta);
    F_inner = local_modal_transform_regls(bright.radius_inner, th, k, m_idx_inner, modalOrder.regEta);
    F_dark  = local_modal_transform_regls(dark.radius,        th, k, m_idx_dark,  modalOrder.regEta);

    spkModel = mr5_proxy_freq(f);

    d_outer_raw = point_source_sum(brightOuterPts, virtual_source_L, virtual_source_R, k);
    targetScale = mean(abs(d_outer_raw));
    d_outer = d_outer_raw / targetScale;

    d_inner_raw = point_source_sum(brightInnerPts, virtual_source_L, virtual_source_R, k);
    d_inner = d_inner_raw / targetScale;

    d_dark = zeros(size(darkBoundaryPts,1), 1);

    d_center_bright = point_source_sum(brightCenterPt, virtual_source_L, virtual_source_R, k) / targetScale;
    d_center_dark = 0;

    beta_outer = F_outer * d_outer;
    beta_inner = F_inner * d_inner;
    beta_dark  = F_dark  * d_dark;

    beta_total = [
        w_outer * beta_outer;
        w_inner * beta_inner;
        darkPenalty * beta_dark
    ];

    B_total = zeros(length(beta_total), Nc);

    for l = 1:Nc
        h_outer_l = reflected_directional_transfer(brightOuterPts, candidates(l), k, spkModel, room);
        h_inner_l = reflected_directional_transfer(brightInnerPts, candidates(l), k, spkModel, room);
        h_dark_l  = reflected_directional_transfer(darkBoundaryPts, candidates(l), k, spkModel, room);

        B_outer_l = F_outer * h_outer_l;
        B_inner_l = F_inner * h_inner_l;
        B_dark_l  = F_dark  * h_dark_l;

        B_total(:,l) = [
            w_outer * B_outer_l;
            w_inner * B_inner_l;
            darkPenalty * B_dark_l
        ];
    end

    beta_total_all{kk} = beta_total;
    B_total_all{kk} = B_total;

    targetData{kk}.d_outer = d_outer;
    targetData{kk}.d_inner = d_inner;
    targetData{kk}.d_dark = d_dark;
    targetData{kk}.d_center_bright = d_center_bright;
    targetData{kk}.d_center_dark = d_center_dark;
    targetData{kk}.targetScale = targetScale;

    modalInfo{kk}.F_outer = F_outer;
    modalInfo{kk}.F_inner = F_inner;
    modalInfo{kk}.F_dark  = F_dark;
    modalInfo{kk}.m_idx_outer = m_idx_outer;
    modalInfo{kk}.m_idx_inner = m_idx_inner;
    modalInfo{kk}.m_idx_dark  = m_idx_dark;
    modalInfo{kk}.beta_outer = beta_outer;
    modalInfo{kk}.beta_inner = beta_inner;
    modalInfo{kk}.beta_dark  = beta_dark;
    modalInfo{kk}.spkModel = spkModel;
end

fprintf('=== 주파수별 local modal block 생성 완료 ===\n');

%% [11] 전 대역 공통 스피커 선택: 순수 local modal objective
fprintf('=== 전 대역 공통 스피커 selection 시작: modal-domain ===\n');

selected = [];
remaining = 1:Nc;
selectionCost = zeros(nSelect,1);

for it = 1:nSelect
    bestJ = inf;
    bestIdx = -1;

    for cand = remaining
        idxTrial = [selected, cand];
        Jsum = 0;

        for kk = 1:K
            Btrial = B_total_all{kk}(:, idxTrial);
            beta_k = beta_total_all{kk};

            qTrial_k = (Btrial' * Btrial + lambda * eye(length(idxTrial))) \ (Btrial' * beta_k);

            % selection 단계는 final 단계의 출력 제약을 강하게 걸지 않음.
            % 후보 위치 전체를 탐색해야 하므로 느슨한 균일 제약만 적용.
            qSafe_trial = qSafe_selection_default * ones(length(idxTrial), 1);
            qTrial_k = apply_weighted_constraints(qTrial_k, qSafe_trial);

            res_k = Btrial * qTrial_k - beta_k;
            Jk = real(res_k' * res_k) + lambda * real(qTrial_k' * qTrial_k);
            Jsum = Jsum + freqWeight(kk) * Jk;
        end

        if Jsum < bestJ
            bestJ = Jsum;
            bestIdx = cand;
        end
    end

    selected = [selected, bestIdx]; %#ok<AGROW>
    remaining(remaining == bestIdx) = [];
    selectionCost(it) = bestJ;

    fprintf('Selection step %d/%d: candidate #%d selected, cost = %.6e\n', ...
        it, nSelect, bestIdx, bestJ);
end

fprintf('=== 스피커 selection 완료 ===\n');

%% [12] 주파수별 final q 계산: 순수 modal-domain + uniform CMP constraint
q_all = cell(K,1);
weightedEnergyUsage = zeros(K,1);
rawArrayEffort = zeros(K,1);

for kk = 1:K
    B_sel_k = B_total_all{kk}(:, selected);
    beta_k = beta_total_all{kk};

    q_k = (B_sel_k' * B_sel_k + lambda * eye(nSelect)) \ (B_sel_k' * beta_k);

    qSafe_selected = qSafe_effective(:, kk);
    q_k = apply_weighted_constraints(q_k, qSafe_selected);

    q_all{kk} = q_k;
    rawArrayEffort(kk) = real(q_k' * q_k);
    weightedEnergyUsage(kk) = sum((abs(q_k).^2) ./ (qSafe_selected.^2));
end

%% [13] 주파수별 성능 평가
outerPressureErr_dB = zeros(K,1);
innerPressureErr_dB = zeros(K,1);
darkEnergy_dB = zeros(K,1);
contrast_dB = zeros(K,1);
centerBrightErr_dB = zeros(K,1);
centerDarkLevel_dB = zeros(K,1);

modalNMSE_outer_dB = zeros(K,1);
modalNMSE_inner_dB = zeros(K,1);
modalResidual_dark_dB = zeros(K,1);
outerCorr = zeros(K,1);
innerCorr = zeros(K,1);

repFreqList = [100 500 1000 1500 2000];
repIdx = zeros(numel(repFreqList),1);
for i = 1:numel(repFreqList)
    [~, repIdx(i)] = min(abs(fList - repFreqList(i)));
end
repTables = struct();

fprintf('=== 주파수별 평가 시작 ===\n');

for kk = 1:K
    f = fList(kk);
    k = 2*pi*f/c;
    q = q_all{kk};
    spkModel = modalInfo{kk}.spkModel;

    H_outer = zeros(size(brightOuterPts,1), nSelect);
    H_inner = zeros(size(brightInnerPts,1), nSelect);
    H_dark  = zeros(size(darkBoundaryPts,1), nSelect);

    for i = 1:nSelect
        H_outer(:,i) = reflected_directional_transfer(brightOuterPts, candidates(selected(i)), k, spkModel, room);
        H_inner(:,i) = reflected_directional_transfer(brightInnerPts, candidates(selected(i)), k, spkModel, room);
        H_dark(:,i)  = reflected_directional_transfer(darkBoundaryPts, candidates(selected(i)), k, spkModel, room);
    end

    P_outer = H_outer * q;
    P_inner = H_inner * q;
    P_dark  = H_dark  * q;

    d_outer = targetData{kk}.d_outer;
    d_inner = targetData{kk}.d_inner;
    d_dark  = targetData{kk}.d_dark;

    % pressure-domain 평가는 참고용이다. 최적화 자체는 modal-domain이다.
    outerPressureErr_dB(kk) = 10*log10(mean(abs(P_outer - d_outer).^2) + 1e-15);
    innerPressureErr_dB(kk) = 10*log10(mean(abs(P_inner - d_inner).^2) + 1e-15);
    darkEnergy_dB(kk) = 10*log10(mean(abs(P_dark).^2) + 1e-15);

    brightEnergyMean = 0.5*(mean(abs(P_outer).^2) + mean(abs(P_inner).^2));
    contrast_dB(kk) = 10*log10((brightEnergyMean + 1e-15)/(mean(abs(P_dark).^2) + 1e-15));

    H_center_bright = zeros(1,nSelect);
    H_center_dark = zeros(1,nSelect);
    for i = 1:nSelect
        H_center_bright(:,i) = reflected_directional_transfer(brightCenterPt, candidates(selected(i)), k, spkModel, room);
        H_center_dark(:,i)   = reflected_directional_transfer(darkCenterPt, candidates(selected(i)), k, spkModel, room);
    end
    P_center_bright = H_center_bright * q;
    P_center_dark = H_center_dark * q;

    d_center_bright = targetData{kk}.d_center_bright;
    d_center_dark = targetData{kk}.d_center_dark;
    centerBrightErr_dB(kk) = 20*log10(abs(P_center_bright - d_center_bright) + 1e-15);
    centerDarkLevel_dB(kk) = 20*log10(abs(P_center_dark - d_center_dark) + 1e-15);

    beta_outer_rec = modalInfo{kk}.F_outer * P_outer;
    beta_inner_rec = modalInfo{kk}.F_inner * P_inner;
    beta_dark_rec  = modalInfo{kk}.F_dark  * P_dark;

    beta_outer = modalInfo{kk}.beta_outer;
    beta_inner = modalInfo{kk}.beta_inner;
    beta_dark  = modalInfo{kk}.beta_dark;

    modalNMSE_outer_dB(kk) = 10*log10(norm(beta_outer_rec - beta_outer)^2/(norm(beta_outer)^2 + 1e-15));
    modalNMSE_inner_dB(kk) = 10*log10(norm(beta_inner_rec - beta_inner)^2/(norm(beta_inner)^2 + 1e-15));
    modalResidual_dark_dB(kk) = 10*log10(norm(beta_dark_rec)^2 + 1e-15);

    outerCorr(kk) = compute_mag_corr(d_outer, P_outer);
    innerCorr(kk) = compute_mag_corr(d_inner, P_inner);

    if any(repIdx == kk)
        fieldName = sprintf('f_%d', f);
        repTables.(fieldName).outer = build_modal_table('outer', modalInfo{kk}.m_idx_outer, beta_outer, beta_outer_rec);
        repTables.(fieldName).inner = build_modal_table('inner', modalInfo{kk}.m_idx_inner, beta_inner, beta_inner_rec);
        repTables.(fieldName).dark  = build_modal_table('dark',  modalInfo{kk}.m_idx_dark,  beta_dark,  beta_dark_rec);
    end
end

fprintf('=== 주파수별 평가 완료 ===\n');

%% [14] 요약 표
summaryTbl = table( ...
    fList(:), ...
    numModes_outer, ...
    numModes_inner, ...
    numModes_dark, ...
    modalNMSE_outer_dB, ...
    modalNMSE_inner_dB, ...
    modalResidual_dark_dB, ...
    outerPressureErr_dB, ...
    innerPressureErr_dB, ...
    darkEnergy_dB, ...
    contrast_dB, ...
    centerBrightErr_dB, ...
    centerDarkLevel_dB, ...
    outerCorr, ...
    innerCorr, ...
    rawArrayEffort, ...
    weightedEnergyUsage, ...
    'VariableNames', ...
    {'FreqHz','NumModesOuter','NumModesInner','NumModesDark', ...
     'ModalNMSE_Outer_dB','ModalNMSE_Inner_dB','ModalResidual_Dark_dB', ...
     'OuterPressureErr_dB','InnerPressureErr_dB','DarkEnergy_dB','Contrast_dB', ...
     'CenterBrightErr_dB','CenterDarkLevel_dB','OuterCorr','InnerCorr', ...
     'RawArrayEffort','WeightedEnergyUsage'});

qSafeSummaryTbl = table((1:nSelect).', ...
    min(qSafe_effective,[],2), mean(qSafe_effective,2), max(qSafe_effective,[],2), ...
    'VariableNames', {'OutputChannel','qSafeMin','qSafeMean','qSafeMax'});

modalOrderSummaryTbl = table(fList(:), numModes_outer, numModes_inner, numModes_dark, ...
    'VariableNames', {'FreqHz','NumModesOuter','NumModesInner','NumModesDark'});

%% [15] 결과 출력
fprintf('\n====================================================\n');
fprintf('순수 로컬 모달 + ceil(kR) adaptive modal order + 균일 CMP 출력 제약 결과\n');
fprintf('주파수 범위: %d Hz ~ %d Hz, step = %d Hz\n', fList(1), fList(end), fList(2)-fList(1));
fprintf('균일 qSafe = %.3f, equivalent Pmax = %.3f\n', qSafe_uniform_value, qSafe_equivalent_Pmax);
fprintf('====================================================\n');

fprintf('선택된 스피커 인덱스:\n');
disp(selected);

fprintf('선택된 스피커 좌표:\n');
for i = 1:nSelect
    fprintf('  #%02d : (%.3f, %.3f), angle = %.1f deg\n', ...
        selected(i), candidates(selected(i)).pos(1), candidates(selected(i)).pos(2), ...
        candidates(selected(i)).angleDeg);
end

fprintf('\n================ 균일 qSafe effective 요약 ================\n');
disp(qSafeSummaryTbl);

fprintf('\n================ 주파수별 adaptive modal mode 개수 앞 20행 ================\n');
disp(modalOrderSummaryTbl(1:min(20,height(modalOrderSummaryTbl)), :));

fprintf('\n================ 전체 주파수 요약 표 앞 20행 ================\n');
disp(summaryTbl(1:min(20,height(summaryTbl)), :));

%% [16] 대표 주파수 1500 Hz 시각화
[~, repK] = min(abs(fList - 1500));
f_rep = fList(repK);
k_rep = 2*pi*f_rep/c;
q_rep = q_all{repK};
spkModel_rep = modalInfo{repK}.spkModel;

xv = linspace(0, room.Lx, Nx);
yv = linspace(0, room.Ly, Ny);
[X, Y] = meshgrid(xv, yv);
gridPts = [X(:), Y(:)];

Pgrid = zeros(size(gridPts,1), 1);
for i = 1:nSelect
    Pgrid = Pgrid + q_rep(i) * reflected_directional_transfer(gridPts, candidates(selected(i)), k_rep, spkModel_rep, room);
end

Pmag2D = reshape(abs(Pgrid), Ny, Nx);
Preal2D = reshape(real(Pgrid), Ny, Nx);

Ptarget = point_source_sum(gridPts, virtual_source_L, virtual_source_R, k_rep) / targetData{repK}.targetScale;
Ptarget2D = reshape(abs(Ptarget), Ny, Nx);

%% [17] 시각화: 배치도
figure('Color','w','Name','Layout');
hold on; axis equal; box on;
xlim([0 room.Lx]); ylim([0 room.Ly]);
rectangle('Position',[0 0 room.Lx room.Ly],'EdgeColor','k','LineWidth',1.5);

allPos = cell2mat(arrayfun(@(s) s.pos, candidates, 'UniformOutput', false)');
plot(allPos(:,1), allPos(:,2), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 9);

selPos = cell2mat(arrayfun(@(s) s.pos, candidates(selected), 'UniformOutput', false)');
plot(selPos(:,1), selPos(:,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);

for i = 1:nSelect
    p = candidates(selected(i)).pos;
    dvec = candidates(selected(i)).dir;
    quiver(p(1), p(2), 0.15*dvec(1), 0.15*dvec(2), 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2);
end

plot(brightOuterPts(:,1), brightOuterPts(:,2), 'b.', 'MarkerSize', 10);
plot(brightInnerPts(:,1), brightInnerPts(:,2), 'c.', 'MarkerSize', 10);
plot(darkBoundaryPts(:,1), darkBoundaryPts(:,2), 'm.', 'MarkerSize', 10);
plot(brightCenterPt(1), brightCenterPt(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(darkCenterPt(1), darkCenterPt(2), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);

draw_circle(bright.center, bright.radius_outer, [0 0.5 1], 1.8);
draw_circle(bright.center, bright.radius_inner, [0 1 1], 1.5);
draw_circle(dark.center, dark.radius, [0.8 0 0.8], 1.8);
plot(virtual_source_L(1), virtual_source_L(2), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
plot(virtual_source_R(1), virtual_source_R(2), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
title('Layout: pure local modal selected speakers');
hold off;

%% [18] 시각화: 주파수별 성능
figure('Color','w','Name','Multiband Metrics');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(fList, modalNMSE_outer_dB, 'LineWidth', 1.8); hold on;
plot(fList, modalNMSE_inner_dB, 'LineWidth', 1.8);
plot(fList, modalResidual_dark_dB, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)'); ylabel('dB');
title('Modal-domain metrics');
legend('Outer modal NMSE','Inner modal NMSE','Dark modal residual','Location','best');

nexttile;
plot(fList, outerPressureErr_dB, 'LineWidth', 1.8); hold on;
plot(fList, innerPressureErr_dB, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)'); ylabel('dB');
title('Pressure mismatch for evaluation');
legend('Outer pressure err','Inner pressure err','Location','best');

nexttile;
plot(fList, darkEnergy_dB, 'LineWidth', 1.8); hold on;
plot(fList, contrast_dB, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)'); ylabel('dB');
title('Dark energy / contrast');
legend('Dark energy','Contrast','Location','best');

nexttile;
plot(fList, weightedEnergyUsage, 'LineWidth', 1.8); hold on;
yline(1.0, 'r--', 'Constraint');
grid on;
xlabel('Frequency (Hz)'); ylabel('Usage');
title('Uniform CMP weighted energy usage');

nexttile;
plot(fList, rawArrayEffort, 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)'); ylabel('||q||^2');
title('Raw array effort');

nexttile;
plot(fList, qSafe_effective(1,:), 'LineWidth', 1.8); hold on;
plot(fList, qSafe_effective(2,:), 'LineWidth', 1.8);
plot(fList, qSafe_effective(3,:), 'LineWidth', 1.8);
plot(fList, qSafe_effective(4,:), 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)'); ylabel('qSafe effective');
title('Uniform effective qSafe by channel');
legend('Ch1','Ch2','Ch3','Ch4','Location','best');

%% [19] 시각화: 1500 Hz target vs reconstructed
figure('Color','w','Name','Target vs Reconstruction at 1500 Hz');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

P_ref = mean(abs(Ptarget));
Ptarget_dB = 20*log10(Ptarget2D/(P_ref + 1e-12) + 1e-12);
Precon_dB = 20*log10(Pmag2D/(P_ref + 1e-12) + 1e-12);

nexttile;
imagesc(xv, yv, Ptarget_dB);
axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
hold on;
draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.5);
draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.2);
draw_circle(dark.center, dark.radius, [1 1 1], 1.5);
plot(virtual_source_L(1), virtual_source_L(2), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
plot(virtual_source_R(1), virtual_source_R(2), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
title(sprintf('Target SPL @ %d Hz', f_rep));
hold off;

nexttile;
imagesc(xv, yv, Precon_dB);
axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
hold on;
draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.5);
draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.2);
draw_circle(dark.center, dark.radius, [1 1 1], 1.5);
plot(selPos(:,1), selPos(:,2), 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title(sprintf('Reconstructed SPL @ %d Hz', f_rep));
hold off;

%% [20] 시각화: 1500 Hz boundary match
H_outer_rep = zeros(size(brightOuterPts,1), nSelect);
H_inner_rep = zeros(size(brightInnerPts,1), nSelect);
H_dark_rep = zeros(size(darkBoundaryPts,1), nSelect);
for i = 1:nSelect
    H_outer_rep(:,i) = reflected_directional_transfer(brightOuterPts, candidates(selected(i)), k_rep, spkModel_rep, room);
    H_inner_rep(:,i) = reflected_directional_transfer(brightInnerPts, candidates(selected(i)), k_rep, spkModel_rep, room);
    H_dark_rep(:,i)  = reflected_directional_transfer(darkBoundaryPts, candidates(selected(i)), k_rep, spkModel_rep, room);
end
P_outer_rep = H_outer_rep*q_rep;
P_inner_rep = H_inner_rep*q_rep;
P_dark_rep = H_dark_rep*q_rep;

d_outer_rep = targetData{repK}.d_outer;
d_inner_rep = targetData{repK}.d_inner;

theta_deg = th*180/pi;
figure('Color','w','Name','Boundary Match at 1500 Hz');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(theta_deg, abs(d_outer_rep), 'k-', 'LineWidth', 2); hold on;
plot(theta_deg, abs(P_outer_rep), 'r--', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Magnitude');
title('Outer boundary magnitude'); legend('Target','Reconstructed','Location','best');

nexttile;
plot(theta_deg, unwrap(angle(d_outer_rep)), 'k-', 'LineWidth', 2); hold on;
plot(theta_deg, unwrap(angle(P_outer_rep)), 'r--', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Phase (rad)');
title('Outer boundary phase'); legend('Target','Reconstructed','Location','best');

nexttile;
plot(theta_deg, abs(d_inner_rep), 'k-', 'LineWidth', 2); hold on;
plot(theta_deg, abs(P_inner_rep), 'r--', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Magnitude');
title('Inner boundary magnitude'); legend('Target','Reconstructed','Location','best');

nexttile;
plot(theta_deg, unwrap(angle(d_inner_rep)), 'k-', 'LineWidth', 2); hold on;
plot(theta_deg, unwrap(angle(P_inner_rep)), 'r--', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Phase (rad)');
title('Inner boundary phase'); legend('Target','Reconstructed','Location','best');


%% [20.5] 추가 시각화: 1500 Hz wrapped phase error
% 기존 Boundary Match figure는 그대로 유지한다.
% 아래 figure는 unwrap branch 차이 때문에 phase가 과도하게 벌어져 보이는지 확인하기 위한 보조 그래프다.
% angle(exp(1i*(phi_rec - phi_target)))은 위상 오차를 [-pi, pi] 범위로 접어서 보여준다.
phaseErr_outer_rep = angle(exp(1i * (angle(P_outer_rep) - angle(d_outer_rep))));
phaseErr_inner_rep = angle(exp(1i * (angle(P_inner_rep) - angle(d_inner_rep))));

phaseErrOuterRMS_rad = sqrt(mean(phaseErr_outer_rep.^2));
phaseErrInnerRMS_rad = sqrt(mean(phaseErr_inner_rep.^2));
phaseErrOuterMax_rad = max(abs(phaseErr_outer_rep));
phaseErrInnerMax_rad = max(abs(phaseErr_inner_rep));

phaseErrSummaryTbl = table( ...
    ["Outer"; "Inner"], ...
    [phaseErrOuterRMS_rad; phaseErrInnerRMS_rad], ...
    [phaseErrOuterMax_rad; phaseErrInnerMax_rad], ...
    'VariableNames', {'Boundary','WrappedPhaseRMSError_rad','WrappedPhaseMaxError_rad'});

fprintf('\n================ 1500 Hz wrapped phase error summary ================\n');
disp(phaseErrSummaryTbl);

figure('Color','w','Name','Wrapped Phase Error at 1500 Hz');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(theta_deg, phaseErr_outer_rep, 'b-', 'LineWidth', 1.8); hold on;
yline(0, 'k--', 'LineWidth', 1.0);
yline(pi, 'r:', 'LineWidth', 1.0);
yline(-pi, 'r:', 'LineWidth', 1.0);
grid on; xlabel('Angle (deg)'); ylabel('Phase error (rad)');
ylim([-pi pi]);
title(sprintf('Outer boundary wrapped phase error @ %d Hz', f_rep));
legend('angle(P_{rec}/P_{target})','0','\pm\pi','Location','best');
hold off;

nexttile;
plot(theta_deg, phaseErr_inner_rep, 'b-', 'LineWidth', 1.8); hold on;
yline(0, 'k--', 'LineWidth', 1.0);
yline(pi, 'r:', 'LineWidth', 1.0);
yline(-pi, 'r:', 'LineWidth', 1.0);
grid on; xlabel('Angle (deg)'); ylabel('Phase error (rad)');
ylim([-pi pi]);
title(sprintf('Inner boundary wrapped phase error @ %d Hz', f_rep));
legend('angle(P_{rec}/P_{target})','0','\pm\pi','Location','best');
hold off;


%% [20.6] 추가 시각화: 1500 Hz gain-scaled boundary match
% 목적:
% - 현재 reconstruction이 target보다 전체적으로 작게 나온 것인지 확인한다.
% - 실제 재생에서 전체 볼륨을 올리는 상황에 해당하도록, 하나의 양의 실수 gain만 사용한다.
% - gain은 outer+inner bright boundary를 함께 기준으로 최소자승 의미에서 구한다.
%
%   g = argmin_{g >= 0} || g * P_bright_rec - P_bright_target ||^2
%
% 이 gain을 outer/inner/dark/전체 grid에 동일하게 적용한다.
% 따라서 magnitude match는 개선될 수 있지만, bright/dark contrast는 원칙적으로 거의 변하지 않는다.
brightTarget_rep = [d_outer_rep; d_inner_rep];
brightRecon_rep  = [P_outer_rep; P_inner_rep];

gainDen = real(brightRecon_rep' * brightRecon_rep) + 1e-15;
g_bright_real = real(brightRecon_rep' * brightTarget_rep) / gainDen;
g_bright_real = max(g_bright_real, 0);

% 참고용 complex scale: 전체 위상 오프셋까지 허용했을 때의 최적 scale
% 실제 볼륨 조절 해석에는 g_bright_real을 사용한다.
alpha_bright_complex = (brightRecon_rep' * brightTarget_rep) / (brightRecon_rep' * brightRecon_rep + 1e-15);

P_outer_gain_rep = g_bright_real * P_outer_rep;
P_inner_gain_rep = g_bright_real * P_inner_rep;
P_dark_gain_rep  = g_bright_real * P_dark_rep;
Pgrid_gain       = g_bright_real * Pgrid;

outerRawNMSE_dB = 10*log10(norm(P_outer_rep - d_outer_rep)^2 / (norm(d_outer_rep)^2 + 1e-15));
outerGainNMSE_dB = 10*log10(norm(P_outer_gain_rep - d_outer_rep)^2 / (norm(d_outer_rep)^2 + 1e-15));
innerRawNMSE_dB = 10*log10(norm(P_inner_rep - d_inner_rep)^2 / (norm(d_inner_rep)^2 + 1e-15));
innerGainNMSE_dB = 10*log10(norm(P_inner_gain_rep - d_inner_rep)^2 / (norm(d_inner_rep)^2 + 1e-15));

brightRawEnergy = 0.5*(mean(abs(P_outer_rep).^2) + mean(abs(P_inner_rep).^2));
brightGainEnergy = 0.5*(mean(abs(P_outer_gain_rep).^2) + mean(abs(P_inner_gain_rep).^2));
darkRawEnergy = mean(abs(P_dark_rep).^2);
darkGainEnergy = mean(abs(P_dark_gain_rep).^2);

contrastRaw_dB = 10*log10((brightRawEnergy + 1e-15)/(darkRawEnergy + 1e-15));
contrastGain_dB = 10*log10((brightGainEnergy + 1e-15)/(darkGainEnergy + 1e-15));
darkRawEnergy_dB = 10*log10(darkRawEnergy + 1e-15);
darkGainEnergy_dB = 10*log10(darkGainEnergy + 1e-15);

% target 대비 평균 magnitude ratio도 같이 출력한다.
outerMeanMagRatio_raw = mean(abs(P_outer_rep)) / (mean(abs(d_outer_rep)) + 1e-15);
outerMeanMagRatio_gain = mean(abs(P_outer_gain_rep)) / (mean(abs(d_outer_rep)) + 1e-15);
innerMeanMagRatio_raw = mean(abs(P_inner_rep)) / (mean(abs(d_inner_rep)) + 1e-15);
innerMeanMagRatio_gain = mean(abs(P_inner_gain_rep)) / (mean(abs(d_inner_rep)) + 1e-15);

gainInfoTbl = table( ...
    g_bright_real, ...
    abs(alpha_bright_complex), ...
    angle(alpha_bright_complex), ...
    'VariableNames', {'RealPositiveGain_BrightBoundary','ComplexGainMagnitude','ComplexGainPhase_rad'});

gainScaleSummaryTbl = table( ...
    ["Outer modal boundary NMSE"; "Inner modal boundary NMSE"; "Dark energy"; "Bright/Dark contrast"; ...
     "Outer mean magnitude ratio"; "Inner mean magnitude ratio"], ...
    [outerRawNMSE_dB; innerRawNMSE_dB; darkRawEnergy_dB; contrastRaw_dB; outerMeanMagRatio_raw; innerMeanMagRatio_raw], ...
    [outerGainNMSE_dB; innerGainNMSE_dB; darkGainEnergy_dB; contrastGain_dB; outerMeanMagRatio_gain; innerMeanMagRatio_gain], ...
    'VariableNames', {'Metric','Raw','GainScaled'});

fprintf('\n================ 1500 Hz gain scaling diagnostic ================\n');
disp(gainInfoTbl);
disp(gainScaleSummaryTbl);

figure('Color','w','Name','Gain-Scaled Boundary Match at 1500 Hz');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(theta_deg, abs(d_outer_rep), 'k-', 'LineWidth', 2); hold on;
plot(theta_deg, abs(P_outer_rep), 'r--', 'LineWidth', 1.6);
plot(theta_deg, abs(P_outer_gain_rep), 'b-.', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Magnitude');
title('Outer boundary magnitude: raw vs gain-scaled');
legend('Target','Raw recon','Gain-scaled recon','Location','best');
hold off;

nexttile;
plot(theta_deg, abs(d_inner_rep), 'k-', 'LineWidth', 2); hold on;
plot(theta_deg, abs(P_inner_rep), 'r--', 'LineWidth', 1.6);
plot(theta_deg, abs(P_inner_gain_rep), 'b-.', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Magnitude');
title('Inner boundary magnitude: raw vs gain-scaled');
legend('Target','Raw recon','Gain-scaled recon','Location','best');
hold off;

nexttile;
plot(theta_deg, abs(P_dark_rep), 'r--', 'LineWidth', 1.6); hold on;
plot(theta_deg, abs(P_dark_gain_rep), 'b-.', 'LineWidth', 1.8);
grid on; xlabel('Angle (deg)'); ylabel('Magnitude');
title('Dark boundary magnitude after same gain');
legend('Raw dark','Gain-scaled dark','Location','best');
hold off;

nexttile;
bar(categorical({'Raw','Gain-scaled'}), [contrastRaw_dB, contrastGain_dB]);
grid on; ylabel('Contrast (dB)');
title('Bright/Dark contrast under global gain');

% 공간 SPL에서도 raw reconstruction과 gain-scaled reconstruction을 같이 비교한다.
PmagGain2D = reshape(abs(Pgrid_gain), Ny, Nx);
PreconGain_dB = 20*log10(PmagGain2D/(P_ref + 1e-12) + 1e-12);

figure('Color','w','Name','Target vs Raw vs Gain-Scaled SPL at 1500 Hz');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(xv, yv, Ptarget_dB);
axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
hold on;
draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.5);
draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.2);
draw_circle(dark.center, dark.radius, [1 1 1], 1.5);
plot(virtual_source_L(1), virtual_source_L(2), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
plot(virtual_source_R(1), virtual_source_R(2), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
title(sprintf('Target SPL @ %d Hz', f_rep));
hold off;

nexttile;
imagesc(xv, yv, Precon_dB);
axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
hold on;
draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.5);
draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.2);
draw_circle(dark.center, dark.radius, [1 1 1], 1.5);
plot(selPos(:,1), selPos(:,2), 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title('Raw reconstruction');
hold off;

nexttile;
imagesc(xv, yv, PreconGain_dB);
axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
hold on;
draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.5);
draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.2);
draw_circle(dark.center, dark.radius, [1 1 1], 1.5);
plot(selPos(:,1), selPos(:,2), 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title(sprintf('Gain-scaled reconstruction, g = %.2f', g_bright_real));
hold off;



%% [20.7] 추가 시각화: 대표 주파수별 결과 비교 (100/500/1000/1500/2000 Hz)
% 기존 1500 Hz 전용 figure는 그대로 유지하고,
% 아래에서는 대표 주파수 100, 500, 1000, 1500, 2000 Hz의 결과를 한 번에 비교한다.
% 목적:
%   1) 저주파/중주파/고주파에서 target vs reconstruction 경향 확인
%   2) 각 주파수에서 wrapped phase error 확인
%   3) global gain scaling이 단순 크기 부족 문제를 얼마나 보정하는지 확인

repFreqCompareList = [100 500 1000 1500 2000];
nRepCompare = numel(repFreqCompareList);

repFreqRows = [];
repFreqPhaseRows = [];

figure('Color','w','Name','Representative Frequencies: Target / Raw / Gain-scaled SPL');
tiledlayout(nRepCompare, 3, 'Padding','compact', 'TileSpacing','compact');

figure('Color','w','Name','Representative Frequencies: Boundary Magnitude');
tiledlayout(nRepCompare, 3, 'Padding','compact', 'TileSpacing','compact');

figure('Color','w','Name','Representative Frequencies: Wrapped Phase Error');
tiledlayout(nRepCompare, 2, 'Padding','compact', 'TileSpacing','compact');

% inner 영역이 boundary인지 pressure point인지 자동 판별
if exist('brightInnerPressurePts','var')
    innerPtsForRep = brightInnerPressurePts;
    innerXBase = 1:size(brightInnerPressurePts,1);
    innerXLabel = 'Inner pressure point index';
    innerTargetStyle = 'ko-';
    innerRawStyle = 'rs--';
    innerGainStyle = 'bd-.';
else
    innerPtsForRep = brightInnerPts;
    innerXBase = theta_deg;
    innerXLabel = 'Angle (deg)';
    innerTargetStyle = 'k-';
    innerRawStyle = 'r--';
    innerGainStyle = 'b-.';
end

for rr = 1:nRepCompare
    f_now = repFreqCompareList(rr);
    [~, kk_now] = min(abs(fList - f_now));
    f_now = fList(kk_now);
    k_now = 2*pi*f_now/c;
    q_now = q_all{kk_now};
    spkModel_now = modalInfo{kk_now}.spkModel;

    H_outer_now = zeros(size(brightOuterPts,1), nSelect);
    H_inner_now = zeros(size(innerPtsForRep,1), nSelect);
    H_dark_now  = zeros(size(darkBoundaryPts,1), nSelect);

    for ii = 1:nSelect
        H_outer_now(:,ii) = reflected_directional_transfer(brightOuterPts, candidates(selected(ii)), k_now, spkModel_now, room);
        H_inner_now(:,ii) = reflected_directional_transfer(innerPtsForRep, candidates(selected(ii)), k_now, spkModel_now, room);
        H_dark_now(:,ii)  = reflected_directional_transfer(darkBoundaryPts, candidates(selected(ii)), k_now, spkModel_now, room);
    end

    P_outer_now = H_outer_now * q_now;
    P_inner_now = H_inner_now * q_now;
    P_dark_now  = H_dark_now  * q_now;

    d_outer_now = targetData{kk_now}.d_outer;
    d_inner_now = targetData{kk_now}.d_inner;

    brightTarget_now = [d_outer_now; d_inner_now];
    brightRecon_now  = [P_outer_now; P_inner_now];
    gainDen_now = real(brightRecon_now' * brightRecon_now) + 1e-15;
    g_now = real(brightRecon_now' * brightTarget_now) / gainDen_now;
    g_now = max(g_now, 0);

    P_outer_gain_now = g_now * P_outer_now;
    P_inner_gain_now = g_now * P_inner_now;
    P_dark_gain_now  = g_now * P_dark_now;

    Pgrid_now = zeros(size(gridPts,1), 1);
    for ii = 1:nSelect
        Pgrid_now = Pgrid_now + q_now(ii) * reflected_directional_transfer(gridPts, candidates(selected(ii)), k_now, spkModel_now, room);
    end
    Pgrid_gain_now = g_now * Pgrid_now;

    Ptarget_now = point_source_sum(gridPts, virtual_source_L, virtual_source_R, k_now) / targetData{kk_now}.targetScale;
    Ptarget2D_now = reshape(abs(Ptarget_now), Ny, Nx);
    Praw2D_now = reshape(abs(Pgrid_now), Ny, Nx);
    Pgain2D_now = reshape(abs(Pgrid_gain_now), Ny, Nx);

    P_ref_now = mean(abs(Ptarget_now));
    Ptarget_dB_now = 20*log10(Ptarget2D_now/(P_ref_now + 1e-12) + 1e-12);
    Praw_dB_now    = 20*log10(Praw2D_now/(P_ref_now + 1e-12) + 1e-12);
    Pgain_dB_now   = 20*log10(Pgain2D_now/(P_ref_now + 1e-12) + 1e-12);

    figure(findobj('Name','Representative Frequencies: Target / Raw / Gain-scaled SPL'));
    nexttile;
    imagesc(xv, yv, Ptarget_dB_now);
    axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
    hold on;
    draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.2);
    draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.0);
    draw_circle(dark.center, dark.radius, [1 1 1], 1.2);
    title(sprintf('Target @ %d Hz', f_now));
    hold off;

    nexttile;
    imagesc(xv, yv, Praw_dB_now);
    axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
    hold on;
    draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.2);
    draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.0);
    draw_circle(dark.center, dark.radius, [1 1 1], 1.2);
    plot(selPos(:,1), selPos(:,2), 'wo', 'MarkerFaceColor','k', 'MarkerSize', 4);
    title(sprintf('Raw recon @ %d Hz', f_now));
    hold off;

    nexttile;
    imagesc(xv, yv, Pgain_dB_now);
    axis xy equal tight; colorbar; colormap turbo; caxis([-25 5]);
    hold on;
    draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.2);
    draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.0);
    draw_circle(dark.center, dark.radius, [1 1 1], 1.2);
    plot(selPos(:,1), selPos(:,2), 'wo', 'MarkerFaceColor','k', 'MarkerSize', 4);
    title(sprintf('Gain-scaled @ %d Hz, g=%.2f', f_now, g_now));
    hold off;

    figure(findobj('Name','Representative Frequencies: Boundary Magnitude'));
    nexttile;
    plot(theta_deg, abs(d_outer_now), 'k-', 'LineWidth', 1.8); hold on;
    plot(theta_deg, abs(P_outer_now), 'r--', 'LineWidth', 1.4);
    plot(theta_deg, abs(P_outer_gain_now), 'b-.', 'LineWidth', 1.5);
    grid on; xlabel('Angle (deg)'); ylabel('|P|');
    title(sprintf('Outer mag @ %d Hz', f_now));
    legend('Target','Raw','Gain','Location','best');
    hold off;

    nexttile;
    plot(innerXBase, abs(d_inner_now), innerTargetStyle, 'LineWidth', 1.8); hold on;
    plot(innerXBase, abs(P_inner_now), innerRawStyle, 'LineWidth', 1.4);
    plot(innerXBase, abs(P_inner_gain_now), innerGainStyle, 'LineWidth', 1.5);
    grid on; xlabel(innerXLabel); ylabel('|P|');
    title(sprintf('Inner mag @ %d Hz', f_now));
    legend('Target','Raw','Gain','Location','best');
    hold off;

    nexttile;
    plot(theta_deg, abs(P_dark_now), 'r--', 'LineWidth', 1.4); hold on;
    plot(theta_deg, abs(P_dark_gain_now), 'b-.', 'LineWidth', 1.5);
    grid on; xlabel('Angle (deg)'); ylabel('|P|');
    title(sprintf('Dark mag @ %d Hz', f_now));
    legend('Raw dark','Gain dark','Location','best');
    hold off;

    phaseErr_outer_now = angle(exp(1i*(angle(P_outer_now) - angle(d_outer_now))));
    phaseErr_inner_now = angle(exp(1i*(angle(P_inner_now) - angle(d_inner_now))));

    phaseOuterRMS_now = sqrt(mean(phaseErr_outer_now.^2));
    phaseInnerRMS_now = sqrt(mean(phaseErr_inner_now.^2));
    phaseOuterMax_now = max(abs(phaseErr_outer_now));
    phaseInnerMax_now = max(abs(phaseErr_inner_now));

    figure(findobj('Name','Representative Frequencies: Wrapped Phase Error'));
    nexttile;
    plot(theta_deg, phaseErr_outer_now, 'b-', 'LineWidth', 1.5); hold on;
    yline(0, 'k--'); yline(pi, 'r:'); yline(-pi, 'r:');
    grid on; ylim([-pi pi]); xlabel('Angle (deg)'); ylabel('rad');
    title(sprintf('Outer phase err @ %d Hz', f_now));
    hold off;

    nexttile;
    plot(innerXBase, phaseErr_inner_now, 'b-', 'LineWidth', 1.5); hold on;
    yline(0, 'k--'); yline(pi, 'r:'); yline(-pi, 'r:');
    grid on; ylim([-pi pi]); xlabel(innerXLabel); ylabel('rad');
    title(sprintf('Inner phase err @ %d Hz', f_now));
    hold off;

    outerRawNMSE_now = 10*log10(norm(P_outer_now - d_outer_now)^2/(norm(d_outer_now)^2 + 1e-15));
    outerGainNMSE_now = 10*log10(norm(P_outer_gain_now - d_outer_now)^2/(norm(d_outer_now)^2 + 1e-15));
    innerRawNMSE_now = 10*log10(norm(P_inner_now - d_inner_now)^2/(norm(d_inner_now)^2 + 1e-15));
    innerGainNMSE_now = 10*log10(norm(P_inner_gain_now - d_inner_now)^2/(norm(d_inner_now)^2 + 1e-15));

    brightRawEnergy_now = 0.5*(mean(abs(P_outer_now).^2) + mean(abs(P_inner_now).^2));
    brightGainEnergy_now = 0.5*(mean(abs(P_outer_gain_now).^2) + mean(abs(P_inner_gain_now).^2));
    darkRawEnergy_now = mean(abs(P_dark_now).^2);
    darkGainEnergy_now = mean(abs(P_dark_gain_now).^2);
    contrastRaw_now = 10*log10((brightRawEnergy_now + 1e-15)/(darkRawEnergy_now + 1e-15));
    contrastGain_now = 10*log10((brightGainEnergy_now + 1e-15)/(darkGainEnergy_now + 1e-15));
    darkRaw_dB_now = 10*log10(darkRawEnergy_now + 1e-15);
    darkGain_dB_now = 10*log10(darkGainEnergy_now + 1e-15);

    repFreqRows = [repFreqRows; ... %#ok<AGROW>
        f_now, g_now, outerRawNMSE_now, outerGainNMSE_now, innerRawNMSE_now, innerGainNMSE_now, ...
        darkRaw_dB_now, darkGain_dB_now, contrastRaw_now, contrastGain_now];

    repFreqPhaseRows = [repFreqPhaseRows; ... %#ok<AGROW>
        f_now, phaseOuterRMS_now, phaseOuterMax_now, phaseInnerRMS_now, phaseInnerMax_now];
end

repFreqGainSummaryTbl = array2table(repFreqRows, ...
    'VariableNames', {'FreqHz','Gain','OuterRawNMSE_dB','OuterGainNMSE_dB', ...
    'InnerRawNMSE_dB','InnerGainNMSE_dB','DarkRawEnergy_dB','DarkGainEnergy_dB', ...
    'ContrastRaw_dB','ContrastGain_dB'});

repFreqPhaseSummaryTbl = array2table(repFreqPhaseRows, ...
    'VariableNames', {'FreqHz','OuterPhaseRMS_rad','OuterPhaseMax_rad', ...
    'InnerPhaseRMS_rad','InnerPhaseMax_rad'});

fprintf('\n================ 대표 주파수 gain scaling summary ================\n');
disp(repFreqGainSummaryTbl);

fprintf('\n================ 대표 주파수 wrapped phase error summary ================\n');
disp(repFreqPhaseSummaryTbl);

%% [21] 시각화: 1500 Hz 실수부
figure('Color','w','Name','Real Part at 1500 Hz');
imagesc(xv, yv, Preal2D);
axis xy equal tight; colorbar; colormap parula; caxis([-2 2]);
hold on;
draw_circle(bright.center, bright.radius_outer, [1 1 1], 1.5);
draw_circle(bright.center, bright.radius_inner, [1 1 1], 1.2);
draw_circle(dark.center, dark.radius, [1 1 1], 1.5);
plot(selPos(:,1), selPos(:,2), 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title(sprintf('Real part Re{P(x,y)} @ %d Hz', f_rep));
hold off;

%% [22] 결과 저장
save('pure_local_modal_adaptive_ceil_uniform_cmp_phasecheck_gaincheck_multirep_2000_result.mat', ...
    'selected', 'q_all', 'summaryTbl', 'qSafe_effective', 'qSafeSummaryTbl', ...
    'qSafe_uniform_value', 'qSafe_equivalent_Pmax', 'modalOrder', ...
    'candidates', 'fList', 'bright', 'dark', 'targetData', 'modalInfo', ...
    'numModes_outer', 'numModes_inner', 'numModes_dark', 'modalOrderSummaryTbl', ...
    'phaseErrSummaryTbl', 'phaseErr_outer_rep', 'phaseErr_inner_rep', ...
    'gainInfoTbl', 'gainScaleSummaryTbl', 'g_bright_real', 'alpha_bright_complex', ...
    'repFreqCompareList', 'repFreqGainSummaryTbl', 'repFreqPhaseSummaryTbl', ...
    'P_outer_gain_rep', 'P_inner_gain_rep', 'P_dark_gain_rep');

fprintf('\n결과 저장 완료: pure_local_modal_adaptive_ceil_uniform_cmp_phasecheck_gaincheck_multirep_2000_result.mat\n');

%% ============================================================
%% 로컬 함수
%% ============================================================

function candidates = generate_arc_candidates(center, radius, anglesDeg)
    candidates = struct('pos', {}, 'dir', {}, 'angleDeg', {});
    for i = 1:numel(anglesDeg)
        th = deg2rad(anglesDeg(i));
        pos = center + radius*[cos(th), sin(th)];
        dirv = center - pos;
        dirv = dirv / norm(dirv);

        candidates(i).pos = pos;
        candidates(i).dir = dirv;
        candidates(i).angleDeg = anglesDeg(i);
    end
end

function p = point_source_sum(points, srcL, srcR, k)
    pts = points;
    rL = sqrt((pts(:,1) - srcL(1)).^2 + (pts(:,2) - srcL(2)).^2);
    rR = sqrt((pts(:,1) - srcR(1)).^2 + (pts(:,2) - srcR(2)).^2);
    rL = max(rL, 0.03);
    rR = max(rR, 0.03);
    p = exp(-1i*k*rL)./(4*pi*rL) + exp(-1i*k*rR)./(4*pi*rR);
end

function m_idx = choose_modal_indices(radius, k, Mmax, besselAbsMin)
% ------------------------------------------------------------
% 주파수별 adaptive modal index 선택
%
% 기본 원칙:
%   M(f,R) = min(Mmax, ceil(kR))
%
% 이유:
% - ceil(kR)는 논문식에서 쓰는 "kR을 넘는 최소 정수" 기준이다.
% - 이전의 ceil(kR+1) 방식은 여기에 margin을 하나 더 붙인 것이어서
%   현재 4스피커 실험에서는 불필요하게 고차 mode를 더 포함할 수 있다.
%
% 추가 안정화:
% - |J_m(kR)|가 besselAbsMin보다 작은 mode는 제외한다.
% - 단, 최소한 하나의 mode는 유지한다. 보통 저주파에서는 m=0만 남는다.
% ------------------------------------------------------------
    x = k * radius;
    Mbase = min(Mmax, max(0, ceil(x)));
    cand = (-Mbase:Mbase).';

    Jabs = abs(besselj(cand, x));
    keep = Jabs >= besselAbsMin;

    if ~any(keep)
        % 모두 제외되면 가장 안정적인 mode 하나를 선택한다.
        [~, id] = max(Jabs);
        keep(id) = true;
    end

    m_idx = cand(keep);
end

function F = local_modal_transform_regls(radius, th, k, m_idx, regEta)
% ------------------------------------------------------------
% boundary pressure p -> local modal coefficient a 변환 행렬 F
%
% 기존 직접식:
%   a_m = (1/N) sum_j p(theta_j) exp(-i m theta_j) / J_m(kR)
%
% 문제:
%   J_m(kR)가 작은 주파수에서 계수가 폭증할 수 있음.
%
% 수정:
%   p(theta_j) ≈ sum_m a_m J_m(kR) exp(i m theta_j)
%   를 regularized least squares로 풀어 a = F p 형태로 만든다.
% ------------------------------------------------------------
    th = th(:);
    Nm = length(m_idx);
    Np = length(th);

    A = zeros(Np, Nm);
    for i = 1:Nm
        m = m_idx(i);
        A(:, i) = besselj(m, k*radius) .* exp(1i*m*th);
    end

    F = (A' * A + regEta * eye(Nm)) \ (A');
end

function qproj = apply_weighted_constraints(q, qSafe)
    qproj = q(:);
    qSafe = qSafe(:);

    for l = 1:length(qproj)
        amp = abs(qproj(l));
        if amp > qSafe(l)
            qproj(l) = qproj(l) * (qSafe(l)/amp);
        end
    end

    eWeighted = sum((abs(qproj).^2) ./ (qSafe.^2));
    if eWeighted > 1 + 1e-12
        qproj = qproj / sqrt(eWeighted);
    end
end

function h_total = reflected_directional_transfer(points, speaker, k, spkModel, room)
    h_direct = base_directional_transfer(points, speaker.pos, speaker.dir, k, spkModel);

    imgPos = [
        -speaker.pos(1),              speaker.pos(2);
         2*room.Lx - speaker.pos(1), speaker.pos(2);
         speaker.pos(1),            -speaker.pos(2);
         speaker.pos(1),             2*room.Ly - speaker.pos(2)
    ];

    imgDir = [
        -speaker.dir(1),  speaker.dir(2);
        -speaker.dir(1),  speaker.dir(2);
         speaker.dir(1), -speaker.dir(2);
         speaker.dir(1), -speaker.dir(2)
    ];

    h_ref = zeros(size(points,1),1);
    for i = 1:4
        h_i = base_directional_transfer(points, imgPos(i,:), imgDir(i,:), k, spkModel);
        h_ref = h_ref + room.refl(i)*h_i;
    end

    h_total = h_direct + h_ref;
end

function h = base_directional_transfer(points, src, dirv, k, spkModel)
    vec = points - src;
    r = sqrt(sum(vec.^2, 2));
    r = max(r, 0.03);

    uv = vec ./ r;
    cosTheta = uv * dirv(:);
    cosTheta = max(min(cosTheta, 1), -1);

    front = max(cosTheta, 0);
    rear = max(-cosTheta, 0);

    D = spkModel.onAxisGain .* ...
        (spkModel.rearFloor + ...
         (1 - spkModel.rearFloor).*(front.^spkModel.pFront) + ...
         spkModel.sideLift.*((1 - abs(cosTheta)).^spkModel.pSide) - ...
         spkModel.rearPenalty.*(rear.^spkModel.pRear));

    D = max(D, spkModel.minGain);
    h = D .* exp(-1i*k.*r) ./ (4*pi*r);
end

function spkModel = mr5_proxy_freq(f)
    % MR5 측정 성향을 단순화한 주파수 의존 프록시 모델
    x = min(max((f - 100)/(2000 - 100), 0), 1);

    spkModel.onAxisGain  = 1.0;
    spkModel.rearFloor   = 0.18 - 0.08*x;
    spkModel.pFront      = 0.8 + 1.2*x;
    spkModel.sideLift    = 0.08 + 0.04*x;
    spkModel.pSide       = 1.0 + 0.4*x;
    spkModel.rearPenalty = 0.04 + 0.08*x;
    spkModel.pRear       = 1.0 + 0.4*x;
    spkModel.minGain     = 0.05;
end

function draw_circle(center, radius, color, lw)
    th = linspace(0, 2*pi, 300);
    x = center(1) + radius*cos(th);
    y = center(2) + radius*sin(th);
    plot(x, y, 'Color', color, 'LineWidth', lw);
end

function rho = compute_mag_corr(a, b)
    a = abs(a(:));
    b = abs(b(:));
    if std(a) < 1e-12 || std(b) < 1e-12
        rho = NaN;
    else
        C = corrcoef(a, b);
        rho = C(1,2);
    end
end

function T = build_modal_table(zoneName, m_idx, betaTarget, betaRecon)
    N = length(m_idx);

    targetReal = real(betaTarget);
    reconReal = real(betaRecon);
    realDiff = reconReal - targetReal;

    targetImag = imag(betaTarget);
    reconImag = imag(betaRecon);
    imagDiff = reconImag - targetImag;

    targetMag = abs(betaTarget);
    reconMag = abs(betaRecon);
    absErr = abs(betaRecon - betaTarget);

    zoneCol = repmat(string(zoneName), N, 1);

    T = table(zoneCol, m_idx, targetReal, reconReal, realDiff, ...
        targetImag, reconImag, imagDiff, targetMag, reconMag, absErr, ...
        'VariableNames', {'Zone','ModeIndex','TargetReal','ReconReal','RealDiff', ...
        'TargetImag','ReconImag','ImagDiff','TargetMag','ReconMag','ComplexAbsError'});
end
