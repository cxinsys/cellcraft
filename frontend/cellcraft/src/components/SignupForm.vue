<template>
  <div>
      <form @submit.prevent="submitForm">
            <div class="field input-field">
              <input type="email" placeholder="Email address" class="input" v-model="email">
            </div>

            <div class="form-link" v-if="!isEmailValidation && email">
              <p>Enter in the email format</p>
            </div>

            <div class="field input-field">
              <input type="name" placeholder="Username" class="input" v-model="username">
            </div>

            <div class="field input-field">
              <input type="password" placeholder="Password" class="password" v-model="password">
            </div>

            <div class="form-link">
              <p>At least 8 characters</p>
            </div>

            <div class="field input-field">
              <input type="password" placeholder="Re-enter password" class="pwCheck" v-model="re_password">
            </div>

            <div class="field button-field">
              <button :disabled="!email || !username || !password" type="submit">Join</button>
            </div>
      </form>
  </div>
</template>

<script>
import { registerUser } from '@/api/index'
import { validateEmail } from '@/utils/validation'

export default {
  data () {
    return {
      email: '',
      password: '',
      username: '',
      re_password: ''
    }
  },
  computed: {
    isEmailValidation () {
      return validateEmail(this.email)
    }
  },
  methods: {
    async submitForm () {
      try {
        const userData = {
          email: this.email,
          password: this.password,
          username: this.username
        }
        console.log(userData)
        const response = await registerUser(userData)
        console.log(response)
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
      this.username = ''
      this.re_password = ''
    }
  }
}
</script>

<style>

</style>
