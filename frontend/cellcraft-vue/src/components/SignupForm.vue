<template>
  <div class="signup-layout">
    <header class="signup-header">
      <h1 class="header-text">Sign up</h1>
    </header>
    <form class="signup-form" @submit.prevent="submitForm">
      <div class="signup__field">
        <input
          class="signup__input"
          type="email"
          placeholder="Email"
          v-model="email"
        />
        <div class="error-message" v-if="!isEmailValidation && email">
          Please follow the email format
        </div>
      </div>

      <div class="signup__field">
        <input
          class="signup__input"
          type="name"
          placeholder="Username"
          v-model="username"
        />
      </div>

      <div class="signup__field">
        <input
          class="signup__input"
          type="password"
          placeholder="Password"
          v-model="password"
        />
        <div class="error-message">At least 8 characters</div>
      </div>

      <div class="signup__field">
        <input
          class="signup__input"
          type="password"
          placeholder="Confirm password"
          v-model="re_password"
        />
      </div>

      <div class="signup__signup">
        <button
          class="signup__button"
          :disabled="!email || !username || !password"
        >
          Join
        </button>
      </div>
    </form>
  </div>
</template>

<script>
import { registerUser } from "@/api/index";
import { validateEmail } from "@/utils/validation";
// import Alert from "@/components/alert";

export default {
  // components: { Alert },
  data() {
    return {
      email: "",
      password: "",
      username: "",
      re_password: "",
      modal: false,
    };
  },
  computed: {
    isEmailValidation() {
      return validateEmail(this.email);
    },
  },
  methods: {
    async submitForm() {
      try {
        const userData = {
          email: this.email,
          password: this.password,
          username: this.username,
        };
        console.log(userData);
        this.modal = true;
        const response = await registerUser(userData);
        console.log(response);
        this.$router.push("/login");
      } catch (error) {
        console.error(error.response.data.detail);
      } finally {
        this.initForm();
      }
    },
    initForm() {
      this.email = "";
      this.password = "";
      this.username = "";
      this.re_password = "";
    },
    async successSignup() {
      this.$router.push("/main");
    },
  },
};
</script>

<style scoped>
.signup-layout {
  width: 100%;
  height: 100%;
  border-radius: 0.4rem;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}
.signup-header {
  width: 100%;
  height: 15%;
  display: flex;
  align-items: center;
  justify-content: center;
  margin-bottom: 5%;
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
.signup-form {
  width: 100%;
  height: 100%;
}
.signup__field {
  width: 100%;
  height: 15%;
  position: relative;
  margin-bottom: 5.5%;
}
.signup__input {
  width: 100%;
  height: 100%;
  border: 1px solid #ccc;
  border-radius: 0.4rem;
  padding: 0 1rem;
  box-sizing: border-box;
  font-size: 1rem;
}
.signup__signup {
  width: 100%;
  height: 15%;
}
.signup__button {
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
.error-message {
  width: 100%;
  height: 5%;
  left: 0.5rem;
  top: calc(100% + 0.3rem);
  position: absolute;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 0.75rem;
  line-height: 0.75rem;
}
</style>
