<template>
  <div>
      <form @submit.prevent="submitForm">
            <div class="field input-field">
              <input type="email" placeholder="이메일을 입력해주세요" class="input" v-model="email">
            </div>

            <div class="form-link" v-if="!isEmailValidation && email">
              <p>이메일 형식으로 입력해주세요</p>
            </div>

            <div class="field input-field">
              <input type="name" placeholder="이름을 입력해주세요" class="input" v-model="username">
            </div>

            <div class="field input-field">
              <input type="password" placeholder="비밀번호를 입력해주세요" class="password" v-model="password">
            </div>

            <div class="form-link">
              <p> 비밀번호 8자리 이상으로 설정해주세요</p>
            </div>

            <div class="field input-field">
              <input type="password" placeholder="비밀번호를 다시 입력해주세요" class="pwCheck" v-model="re_password">
            </div>

            <div class="field button-field">
              <button :disabled="!email || !username || !password" type="submit">가입하기</button>
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
