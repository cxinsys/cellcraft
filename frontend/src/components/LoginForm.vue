<template>
  <div class="login-layout">
    <header class="login-header">
      <h1 class="header-text">Login</h1>
    </header>
    <form class="login-form" @submit.prevent="submitForm">
      <div class="login__field">
        <input
          class="login__input"
          type="email"
          placeholder="Email"
          v-model="email"
        />
      </div>

      <div class="login__field">
        <input
          class="login__input"
          type="password"
          placeholder="Password"
          v-model="password"
        />
      </div>

      <div class="login__error" v-if="isError">{{ errorMessage }}</div>

      <div class="login__login">
        <button
          class="login__button"
          :disabled="!email || !password"
          type="submit"
        >
          Login
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
        console.log(userData);
        await this.$store.dispatch("LOGIN", userData);
        this.$router.push("/projects");
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
  display: flex;
  align-items: center;
  justify-content: center;
}
.header-text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 2rem;
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
  height: 25%;
  margin-top: 2.5%;
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
  height: 25%;
  margin-top: 5%;
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
