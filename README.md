# 0️⃣Pre-ready
If you want to use this repository, need to install WSL2 on Windows. 
This recommand Ubuntu24-04 LTS version.
On Linux system(**NOT** WSL), naturaly can run.

# 1️⃣Make binary
Push below command on your CMD(bash), respectively.

``` bash
  cmake .
```

```bash
  make
```

# 2️⃣Run
```bash
  ./main
```

# 3️⃣task(20260117-0832 updated)
## 완료됨
### 260112
1. 벡터, 행렬, 암호문 구조화
2. 연산 형태 구조화
3. 필요한 연산기 구조화
### 260113
1. 구현
### 260114
1. 모든 mpz_t 타입을 c++ mpz_class로 변경
2. class에 담기는 배열을 std::vector로 변경
3. 현재 구현된 원시근 구하는것, 다항식 ntt 연산 부분 정리
4. cipher 클래스 재수정(거의 대부분)
5. 메모리 최적화
6. 암호 합, 곱 연산 추가
7. SEAL 없애버리고 batch_encoder 구현 추가(2^15보다 큰 수도 가능:[느림])
### 260115
1. ARX 제어기 연산 구현
2. 암호문 재활용 방안 구현
3. 지수 연산 구현
### 260117
1. 암호문의 행렬로의 재표현 구현
~2. 행렬-벡터 곱과 같은 형태의 제어기 구현~
3. 암호문의 행렬표현의 메모리 사용량이 너무 커서 EKF만 뽑을때만 쓰고 실제 연산에서 사용하지 않도록 수정(기존 2^14 ring dimension 기준 matrix 표현 하나당 8.5GB | arx 계수 P Q에 대한 암호문 8개 기준 68GB)

## 할것
1. VC 구현
2. 공격 시나리오 설정
