<template>
  <div class="login-layout">
    <header class="login-header">
      <h1 class="header-text">Change Password</h1>
    </header>
    <p class="login__notice" v-if="isRequired">
      You must change your password before continuing.
    </p>
    <form class="login-form" @submit.prevent="submitForm">
      <div class="login__field">
        <input class="login__input" type="password" placeholder="Current password"
          autocomplete="current-password" v-model="currentPassword" />
      </div>
      <div class="login__field">
        <input class="login__input" type="password" :placeholder="`New password (${minLength}+ characters)`"
          autocomplete="new-password" v-model="newPassword" />
      </div>
      <div class="login__field">
        <input class="login__input" type="password" placeholder="Confirm new password"
          autocomplete="new-password" v-model="confirmPassword" />
      </div>
      <div class="login__error" v-if="errorMessage">{{ errorMessage }}</div>
      <div class="login__success" v-if="successMessage">{{ successMessage }}</div>
      <div class="login__login">
        <button
          class="login__button"
          :disabled="!canSubmit"
          type="submit"
          :class="{ 'login__button--loading': isLoading }"
        >
          {{ isLoading ? 'Changing...' : 'Change Password' }}
        </button>
      </div>
    </form>
  </div>
</template>

<script>
import { changePassword } from "@/api/index";

export const MIN_PASSWORD_LENGTH = 8;

export default {
  data() {
    return {
      currentPassword: "",
      newPassword: "",
      confirmPassword: "",
      errorMessage: "",
      successMessage: "",
      isLoading: false,
      minLength: MIN_PASSWORD_LENGTH,
    };
  },
  computed: {
    isRequired() {
      return this.$store.getters.mustChangePassword;
    },
    canSubmit() {
      return this.currentPassword && this.newPassword && this.confirmPassword && !this.isLoading;
    },
  },
  methods: {
    validate() {
      if (this.newPassword.length < MIN_PASSWORD_LENGTH) {
        return `New password must be at least ${MIN_PASSWORD_LENGTH} characters.`;
      }
      if (this.newPassword !== this.confirmPassword) {
        return "New passwords do not match.";
      }
      if (this.newPassword === this.currentPassword) {
        return "New password must be different from the current password.";
      }
      return "";
    },
    async submitForm() {
      this.errorMessage = this.validate();
      this.successMessage = "";
      if (this.errorMessage) return;

      this.isLoading = true;
      try {
        await changePassword(this.currentPassword, this.newPassword);
        this.$store.commit("setPasswordChanged");
        this.successMessage = "Password changed.";
        this.currentPassword = "";
        this.newPassword = "";
        this.confirmPassword = "";
        this.$router.push(this.$store.getters.isSuperUser ? "/admin" : "/projects");
      } catch (error) {
        const detail = error.response && error.response.data && error.response.data.detail;
        this.errorMessage = typeof detail === "string" ? detail : "Failed to change password. Please try again.";
      } finally {
        this.isLoading = false;
      }
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
  padding-bottom: 1.5rem;
}

.login-header {
  width: 90%;
  margin: 6% 0 4% 0;
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
  color: rgb(81, 81, 81);
}

.login-form {
  width: 90%;
}

.login__field {
  width: 100%;
  height: 3rem;
  margin-top: 0.75rem;
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
  height: 3rem;
  margin-top: 1.5rem;
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
  font-weight: 500;
  font-size: 1rem;
}

.login__button:hover:not(:disabled) {
  background: rgb(65, 109, 199);
}

.login__button:disabled {
  background: rgb(170, 170, 170);
  cursor: not-allowed;
  opacity: 0.7;
}

.login__notice,
.login__error,
.login__success {
  width: 90%;
  margin: 0.75rem 0 0 0;
  font-family: "Montserrat", sans-serif;
  font-weight: 500;
  font-size: 0.95rem;
  line-height: 1.3rem;
}

.login__notice {
  color: rgb(81, 81, 81);
}

.login__error {
  color: red;
}

.login__success {
  color: rgb(40, 140, 80);
}
</style>
