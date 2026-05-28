# Local Modal Sound Field Control for Bright/Dark Zone Reproduction

This project implements a MATLAB-based simulation framework for multizone sound field control using local modal domain analysis. The goal is to reproduce a desired sound field in a bright zone while suppressing sound energy in a dark zone using a selected loudspeaker array.

## 1. Project Motivation
개인 공간의 음장 제어는 과거부터 계속 되어온 관심사입니다. 이를 해결하기 위한 시도는 계속 발전되어왔습니다.
그 중 기존의 음압 기반 제어(Pressure-point matching)은 제어점 사이의 빈 공간에서 음장이 심하게 왜곡되는 문제가 발생하며, 방 전체를 커버하는 글로벌 모달(Global Modal) 방식은 주파수가 높아질수록 연산량이 기하급수적으로 폭발하는 한계가 있습니다.

본 프로젝트는 개인 공간이라는 점에 주목하여 이 한계를 극복하기 위해 방전체가 아닌 귀 주변의 좁은 영역(Local)으로 계산 구역을 좁히고, 모달 차수를 물리적 한계인 $N = \lceil kR \rceil$로 제한하여 연산량을 낮추면서도 안정적인 음장을 만들어내는 순수 Local modal 제어 시뮬레이터를 구현했습니다.
나아가 현재 실제 공간에서의 측정 실험을 통해 실제 공간에서의 유효성을 판단하고자 합니다.

## 2. Methodology & Optimization Tactics
본 시뮬레이션 다음 상황을 가정합니다.
  1. Direct-Path Acoustic Propagation modeled by Free-Field Green's Function
  2. First-Order Early Reflection Modeling via Image Source Method (ISM)
  3. Multizone Sound Field Control (SFC) based on Spatial Audio Reproduction
논문 Optimizing loudspeaker placement based on modal domain analysis for multizone sound field reproduction 의 내용을 참고하여 시뮬레이션 코드를 작성하였습니다.
본 시뮬레이션은 결과 대조를 위해 다음 두 가지 모델을 비교 분석합니다.

### Proposed : Pure Local Modal Domain Method
  * Dark zone (소멸): 헬름홀츠 방정식의 경계값 정리를 차용하여 1개의 Boundary만 설정한다. 다만 내부 측정을 위해 중앙에 측정 점 하나를 설정한다.
  * Bright zone (구현): 목표 음장의 Local한 구현을 위해 Inner/Outer 2개의 Boundary를 설정한다. 다만 내부 측정을 위해 중앙에 측정 점 하나를 설정한다.

### Control Group: Local-Pressure Modal Hybrid (4-points)
  * 대조군으로, 위 방식에서 Outer Boundary 내부에 Inner boundary 가 아닌 4개의 음압점을 균등하게 잡아 매칭 타겟으로 삼는다. 이후 local 연산을 수행한다.

## 3. Main Features

- Pure local modal-domain sound field control
- Local-pressure modal hybrid comparison model
- Adaptive modal order based on \(N = \lceil kR \rceil\)
- Uniform CMP-style output constraint
- Wrapped phase error analysis
- Gain-scaled reconstruction diagnostic
- Representative frequency comparison at 100, 500, 1000, 1500, and 2000 Hz

## 4. Results

Representative results are provided in the `results/figures` directory.

The simulation compares target, raw reconstruction, and gain-scaled reconstruction across multiple representative frequencies. Wrapped phase error is also used to avoid misleading interpretation caused by phase unwrapping.

<img width="2559" height="1219" alt="image" src="https://github.com/user-attachments/assets/4eaf3f2a-f02a-41be-bce3-a83d281c3b9d" />

<img width="2552" height="1214" alt="image" src="https://github.com/user-attachments/assets/4ee69693-cd49-4299-ba53-07502901f8a8" />
: pure local modal method result

<img width="2559" height="1221" alt="image" src="https://github.com/user-attachments/assets/a6df65fe-42e7-47c7-ad65-e2e1b987e653" />

<img width="1132" height="1109" alt="image" src="https://github.com/user-attachments/assets/7a775ca8-cb58-4f04-a8d2-f28879381d8e" />
: pressure local modal method result
