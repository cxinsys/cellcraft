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
        <button class="login__button" :disabled="!email || !password" type="submit">
          Sign In
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
    };
  },
  methods: {
    async submitForm() {
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
      } catch (error) {
        console.error(error.response.data.detail);
        this.errorMessage = error.response.data.detail;
        this.isError = true;
      } finally {
        this.initForm();
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

  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  text-decoration: none;
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
