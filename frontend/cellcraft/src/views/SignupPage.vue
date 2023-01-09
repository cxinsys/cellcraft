<template>
  <div class="layout">
    <div class="signup-component">
        <router-link class="signup__logo" to="/main">
          <img class="logo__img" src="@/assets/cellcraft_logo.png">
          <p class="logo__text">CELLCRAFT</p>
        </router-link>
        <div class="signup-box">
          <SignupForm></SignupForm>
        </div>
        <footer class="signup-footer">
          <div class="copyright__txt">
            CELLCRAFT © 2023. All rights reserved
          </div>
        </footer>
    </div>
  </div>
</template>

<script>
import SignupForm from '@/components/SignupForm.vue'

export default {
  components: {
    SignupForm
  }
}
</script>

<style scoped>
 .layout{
  width: 100%;
  height: 100%;
  display: flex;
  align-items: flex-start;
  justify-content: center;
  padding: 5rem;
  box-sizing: border-box;
}
.signup-component{
  width: 26rem;
  height: 40rem;
}
.signup__logo{
  width: 100%;
  height: 6rem;
  display: flex;
  justify-content: center;
  align-items: center;
  text-decoration: none;
  color: black;
  margin-bottom: 1rem;
}
.logo__img{
  width: 6rem;
  height: 6rem;
  object-fit: contain;
}
.logo__text{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 3rem;
  line-height: 3rem;
  text-decoration: none;
}
.signup-box{
  width: 100%;
  height: 30rem;
  margin-bottom: 2rem;
}
.signup-footer{
  width: 100%;
  height: 3rem;
  border-top: 1px solid rgba(204, 204, 204, 0.5);
}
</style>
