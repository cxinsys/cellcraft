<template>
  <!-- <div class="layout">
    <form class="fileUpload-form" @submit.prevent="uploadFile">
        <label class="fileUpload-form__title">File</label>
        <div class="toggle-layout" v-if="toggle_file">
          <div class="toggle-layout__row">
            <label class="fileUpload-form__button">
              파일 찾기
              <input df-file class="fileUpload-form__input" type="file" ref="selectFile" @change.prevent="previewFile" accept="text/csv" />
            </label>
            <label class="fileUpload-form__button">
              업로드
              <input class="fileUpload-form__button" type="submit" value="업로드">
            </label>
          </div>
          <input class="fileUpload-form__info" placeholder="첨부파일" v-if="selectFile" v-model="selectFile.name" readonly>
          <input class="fileUpload-form__info" placeholder="첨부파일" v-else readonly>
          <ul class="fileUploa d-form__info" v-if="selectFile">
            <li>name : {{ selectFile.name }}</li>
            <li>size : {{ selectFile.size | formatBytes}}</li>
            <li>type : {{ selectFile.type }}</li>
          </ul>
          <ul class="fileUpload-form__info" v-else>
            <li>name : </li>
            <li>size : </li>
            <li>type : </li>
          </ul>
        </div>
    </form>
    <div class="toggle-icon" @click="toggleFile">
      <i class="fileUpload-form__arrow fa-solid fa-arrow-up" v-if="toggle_arrow"></i>
      <i class="fileUpload-form__arrow fa-solid fa-arrow-down" v-else></i>
    </div>
  </div> -->
  <div>
    <div class="nodeBox">
      <img class="nodeBox__icon" src="@/assets/file-upload.png">
    </div>
  </div>
</template>

<script>
import { uploadForm } from '@/api/index'

export default {
  data () {
    return {
      title: 'File Upload',
      selectFile: null,
      is_upload: false,
      done_upload: false,
      toggle_arrow: false,
      toggle_file: false
    }
  },
  methods: {
    previewFile () {
      if (this.$refs.selectFile.files.length > 0) {
        this.selectFile = this.$refs.selectFile.files[0]
      }
    },
    toggleFile () {
      this.toggle_arrow = !this.toggle_arrow
      this.toggle_file = !this.toggle_file
    },
    async uploadFile () {
      if (this.selectFile) {
        this.is_upload = true
        const form = new FormData()
        form.append('files', this.selectFile)
        const response = await uploadForm(form)
        console.log(response)
        if (response) {
          this.done_upload = true
        }
      }
    }
  },
  filters: {
    formatBytes (a, b) {
      if (a === 0) return '0 Bytes'
      const c = 1024
      const d = b || 2
      const e = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB']
      const f = Math.floor(Math.log(a) / Math.log(c))

      return parseFloat((a / Math.pow(c, f)).toFixed(d)) + ' ' + e[f]
    }
  }
}
</script>

<style scoped>
.nodeBox{
  width: 3rem;
  height: 3rem;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
}
.nodeBox__icon{
  width: 70%;
  height: 70%;
  object-fit: contain;
  filter: invert(97%) sepia(99%) saturate(0%) hue-rotate(123deg) brightness(107%) contrast(101%);
}
.layout{
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}
.toggle-layout{
  margin: 20px;
  border-radius: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}
.toggle-layout__row{
  width: 100%;
  height: auto;
  display: flex;
  align-items: center;
  justify-content:space-between;
}
.fileUpload-form{
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
}
.fileUpload-form__button{
  display: inline-block;
  width: 46%;
  height: 1.5rem;
  padding-top: 10px;
  border: 1px solid black;
  color: black;
  background: white;
  margin-bottom: 5px;
  cursor: pointer;
  vertical-align: middle;
  font-size: 14px;
}

.fileUpload-form input[type="file"], .fileUpload-form input[type="submit"] {
    position: absolute;
    width: 0;
    height: 0;
    padding: 0;
    overflow: hidden;
    border: 0;
}

.fileUpload-form__title{
    font-size: 1rem;
    font-weight: bold;
}

/* .fileUpload-form__info{

} */

.fileUpload-form__arrow{
  color: black;
  width: 1.5rem;
  height: auto;
}
.toggle-icon{
  background: rgb(227, 243, 252);
  width: 20px;
  height: 20px;
  border-radius: 50%;
  margin-top: 5px;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
}
</style>
