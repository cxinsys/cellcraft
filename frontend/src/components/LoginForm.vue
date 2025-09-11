<template>
  <div class="login-layout">
    <header class="login-header">
      <h1 class="header-text">Sign in to Cellcraft</h1>
    </header>
    <form class="login-form" @submit.prevent="submitForm">
      <div class="login__field">
        <input class="login__input" type="email" placeholder="Email" v-model="email" />
      </div>

      <div class="login__field">
        <input class="login__input" type="password" placeholder="Password" v-model="password" />
      </div>

      <div class="login__error" v-if="isError">{{ errorMessage }}</div>

      <div class="login__login">
        <button 
          class="login__button" 
          :disabled="!email || !password || isLoading" 
          type="submit"
          :class="{ 'login__button--loading': isLoading }"
        >
          {{ isLoading ? 'Signing In...' : 'Sign In' }}
        </button>
      </div>
    </form>
  </div>
</template>

<script>
export default {
  data() {
    return {
      email: "",
      password: "",
      errorMessage: "",
      isError: false,
      isLoading: false,
    };
  },
  methods: {
    async submitForm() {
      // 로딩 시작 및 에러 상태 클리어
      this.isLoading = true;
      this.isError = false;
      this.errorMessage = "";

      try {
        const userData = {
          username: this.email,
          password: this.password,
        };
        await this.$store.dispatch("LOGIN", userData);

        // isSuperUser 여부 확인 - 수정된 로직
        await this.$nextTick(); // 상태 업데이트 대기
        const isSuperUser = this.$store.getters.isSuperUser;
        
        if (isSuperUser === true) {
          this.$router.push("/admin");
        } else {
          this.$router.push("/projects");
        }

        // 성공 시에만 폼 초기화
        this.initForm();
      } catch (error) {
        // 개발자용 콘솔 로그
        console.error("로그인 오류:", error);

        // 사용자에게 표시할 오류 메시지 결정
        let userErrorMessage = "Login failed. Please try again later.";

        if (error.response) {
          // 서버 응답이 있는 경우
          if (error.response.status === 401) {
            userErrorMessage = "Invalid email or password.";
          } else if (error.response.status === 422) {
            userErrorMessage = "Invalid input information.";
          } else if (error.response.status >= 500) {
            userErrorMessage = "Server error occurred. Please try again later.";
          } else if (error.response.data && error.response.data.detail) {
            // 서버에서 제공하는 상세 메시지가 있는 경우
            userErrorMessage = error.response.data.detail;
          }
        } else if (error.request) {
          // 네트워크 오류
          userErrorMessage = "Please check your network connection.";
        }

        this.errorMessage = userErrorMessage;
        this.isError = true;
      } finally {
        // 로딩 상태만 해제 (오류 발생 시 폼 유지)
        this.isLoading = false;
      }
    },
    initForm() {
      this.email = "";
      this.password = "";
    },
  },
};
</script>

<style scoped>
.login-layout {
  width: 100%;
  height: 100%;
  background: rgba(204, 218, 245, 0.6);
  border-radius: 0.4rem;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.login-header {
  width: 90%;
  height: 15%;
  margin: 8% 0 7% 0;
  display: flex;
  align-items: center;
  justify-content: center;
}

.header-text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.3rem;
  line-height: 2rem;
  text-decoration: none;
  color: rgb(81, 81, 81);
}

.login-form {
  width: 90%;
  height: 80%;
}

.login__field {
  width: 100%;
  height: 20%;
  margin-top: 1.5%;
  position: relative;
}

.login__input {
  width: 100%;
  height: 100%;
  border: 1px solid #ccc;
  border-radius: 0.4rem;
  padding: 0 1rem;
  box-sizing: border-box;
  font-size: 1rem;
}

.login__input:focus {
  border: 1px solid rgb(40, 84, 197);
}

.login__login {
  width: 100%;
  height: 15%;
  margin-top: 15%;
}

.login__button {
  width: 100%;
  height: 100%;
  border-radius: 0.4rem;
  background: rgb(75, 119, 209);
  color: white;
  border: none;
  cursor: pointer;
  transition: background-color 0.3s ease, opacity 0.3s ease;

  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  text-decoration: none;
}

.login__button:hover:not(:disabled) {
  background: rgb(65, 109, 199);
}

.login__button:disabled {
  background: rgb(170, 170, 170);
  cursor: not-allowed;
  opacity: 0.7;
}

.login__button--loading {
  position: relative;
}

.login__error {
  width: 100%;
  margin: 1rem 0 0.5rem 0;
  display: flex;
  align-items: center;
  color: red;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  text-decoration: none;
}
</style>
