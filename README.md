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

# 3️⃣task(20260114-0606 updated)
1. 모든 mpz_t 타입을 c++ mpz_class로 변경 (v)
2. class에 담기는 배열을 std::vector로 변경 (v)
3. 현재 구현된 원시근 구하는것, 다항식 ntt 연산 부분 정리 (v)
4. cipher 클래스 재수정(거의 대부분) (v)
5. 메모리 최적화 (v)
6. 암호 합, 곱 연산 추가 (v)

