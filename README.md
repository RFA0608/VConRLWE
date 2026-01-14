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

# 3️⃣task(20260115-0703 updated)
## 완료됨
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

## 할것
3. VC 구현
4. 공격 시나리오 설정

