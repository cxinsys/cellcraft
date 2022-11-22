<style scoped>
.modal__alert__header{
  display: flex;
  align-content: center;
  justify-content: flex-end;
  align-items: center;
  padding: 10px;
}
.modal__alert__header__closeBtn:hover{
  cursor: default;
}
.modal__alert__main{
  display: flex;
  flex-direction: column;
  align-content: center;
  justify-content: center;
  align-items: center;
  margin-top: 20px;
}
.modal__alert__main__content{
  /* margin-top: 50px; */
  position: absolute;
  top: 45%;
  display: flex;
  flex-direction: column;
  align-items: center;
}
.modal__alert__main__btn{
  margin-top: 50px;
}
</style>

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

        <!-- <div class="field button-field" > -->
        <div class="field button-field">
          <button :disabled="!email || !username || !password">Join</button>
        </div>
        <Alert @click="modal=false" v-if="modal">
          <div slot="header" class="modal__alert__header">
          </div>
          <div slot="body" class="modal__alert__main">
            <div class="modal__alert__main__content">
              Success Signup!<br><br>Move to MainPage
              <button class="modal__alert__main__btn" @click="successSignup">OK</button>
            </div>
            <!-- <button class="modal__alert__main__btn" type="submit">OK</button> -->
          </div>
        </Alert>
      </form>

  </div>
</template>

<script>
import { registerUser } from '@/api/index'
import { validateEmail } from '@/utils/validation'
import Alert from '@/components/alert'

export default {
  components: {Alert},
  data () {
    return {
      email: '',
      password: '',
      username: '',
      re_password: '',
      modal: false
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
        this.modal = true
        const response = await registerUser(userData)
        console.log(response)
        // this.$router.push('/main')
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
    },
    async successSignup () {
      this.$router.push('/main')
    }
  }
}
</script>
