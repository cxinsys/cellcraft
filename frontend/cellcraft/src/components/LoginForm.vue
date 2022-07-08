<template>
  <div>
      <form @submit.prevent="submitForm">
            <div class="field input-field">
              <input type="email" placeholder="Email" class="input" v-model="email" >
            </div>

            <div class="field input-field">
              <input type="password" placeholder="Password" class="password" v-model="password" >
            </div>

            <div class="form-link">
              <a href="#" class="forgot-pass">Forgot password?</a>
            </div>

            <div class="field button-field">
              <button :disabled="!email || !password" type="submit">Login</button>
            </div>
      </form>
  </div>
</template>

<script>
export default {
  data () {
    return {
      email: '',
      password: ''
    }
  },
  methods: {
    async submitForm () {
      try {
        const userData = {
          username: this.email,
          password: this.password
        }
        await this.$store.dispatch('LOGIN', userData)
        this.$router.push('/main')
      } catch (error) {
        console.error(error.response.data.detail)
      } finally {
        this.initForm()
      }
    },
    initForm () {
      this.email = ''
      this.password = ''
    }
  }
}
</script>

<style>

</style>
