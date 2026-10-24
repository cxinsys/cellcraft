<template>
  <div class="layout">
    <div class="login-component">
      <router-link class="login__logo" to="/">
        <img class="logo__img" src="@/assets/cellcraft_logo.png" />
        <p class="logo__text">CELLCRAFT</p>
      </router-link>
      <div class="login-box">
        <ChangePasswordForm></ChangePasswordForm>
      </div>
    </div>
  </div>
</template>

<script>
import ChangePasswordForm from "@/components/ChangePasswordForm.vue";

export default {
  components: {
    ChangePasswordForm,
  },
};
</script>

<style scoped>
.layout {
  width: 100vw;
  height: 100vh;
  display: flex;
  align-items: flex-start;
  justify-content: center;
  box-sizing: border-box;
}

.login-component {
  width: 26rem;
  margin-top: 5rem;
}

.login__logo {
  width: 100%;
  height: 6rem;
  display: flex;
  justify-content: center;
  align-items: center;
  text-decoration: none;
  color: black;
  margin-bottom: 2rem;
}

.logo__img {
  width: 6rem;
  height: 6rem;
  object-fit: contain;
}

.logo__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 800;
  font-size: 3rem;
  color: rgb(40, 84, 197);
  line-height: 3rem;
}

.login-box {
  width: 100%;
}
</style>
