<template>
  <div class="layout">
    <form class="fileUpload-form" @submit.prevent="uploadFile">
      <div class="form-row">
        <h1 class="form__name">Choose File</h1>
        <div class="form__selectFile">
          <ul class="form__fileList">
            <li class="fileList__item" v-for="(file, idx) in files_list" :key="idx" @click="fileSelect">
              <p class="fileList__text">{{file.file_name}}</p>
            </li>
            <li class="fileList__item" v-if="files_list.length === 0">
              <p class="fileList__text--blank">Please upload a new file</p>
            </li>
          </ul>
          <label class="form__button">
            <div class="form__addfile">
              <img class="form__button--foldericon" src="@/assets/add-file.png">
              <img class="form__button--plusicon" src="@/assets/plus.png">
            </div>
            <input df-file class="form__input" type="file" ref="selectFile" @change.prevent="previewFile" accept="text/csv" />
          </label>
        </div>
      </div>
      <div class="form-row">
        <h1 class="form__name">Current File</h1>
        <ul class="form__info" v-if="selectFile">
          <li class="form__info--type">
            <img class="form__info--img" src="@/assets/csv.png">
          </li>
          <li class="form__info--name">{{ selectFile.name }}&nbsp;&nbsp;&nbsp;{{ selectFile.size | formatBytes}}</li>
        </ul>
        <ul class="form__info" v-else>
          <li class="form__info--blank">Please add data file</li>
        </ul>
      </div>
      <div class="form-row">
        <label class="form__button--apply">
          Apply
          <input class="form__input" type="submit" value="업로드">
        </label>
      </div>
    </form>
  </div>
</template>

<script>
import { uploadForm, getFiles, findFile } from '@/api/index'

export default {
  data () {
    return {
      node_name: 'File',
      selectFile: null,
      is_upload: false,
      done_upload: false,
      files_list: []
    }
  },
  methods: {
    previewFile () {
      if (this.$refs.selectFile.files.length > 0) {
        this.selectFile = this.$refs.selectFile.files[0]
      }
    },
    async uploadFile () {
      if (this.selectFile === this.$refs.selectFile.files[0]) {
        this.is_upload = true
        const form = new FormData()
        form.append('files', this.selectFile)
        const response = await uploadForm(form)
        console.log(response)
        if (response) {
          this.done_upload = true
          const fileList = await getFiles()
          console.log(fileList.data)
          this.files_list = fileList.data
        }
      }
    },
    async fileSelect (ev) {
      console.dir(ev.target.innerText)
      const fileInfo = {
        file_name: ev.target.innerText
      }
      const selectFile = await findFile(fileInfo)
      this.selectFile = {
        name: selectFile.data.file_name,
        size: selectFile.data.file_size
      }
      console.log(selectFile.data)
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
  },
  async mounted () {
    const fileList = await getFiles()
    console.log(fileList.data)
    this.files_list = fileList.data
  }
}
</script>

<style scoped>
.layout{
  width: 100%;
  height: 100%;
  position: relative;
}
.fileUpload-form{
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
}
.form-row{
  width: 100%;
  height: 30%;
  position: relative;
}
.form-row:first-child{
  height: 60%;
}
.form-row:last-child{
  height: 10%;
}
.form__name{
  width: 100%;
  height: 30%;
  display: flex;
  align-items: center;
  padding: 0 2rem;
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 2rem;
  line-height: 2rem;
  color: rgb(51, 51, 51);
}
.form__selectFile{
  width: 100%;
  height: 70%;
  display: flex;
}
.form__fileList{
  width: 60%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  overflow: hidden;
}
.fileList__item{
  width: 70%;
  padding: 0.5rem;
  background: #FFFFFF;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  margin-bottom: 1rem;
}
.fileList__text{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.3rem;
  line-height: 1.3rem;
  color: rgb(51, 51, 51);
}
.fileList__text--blank{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.3rem;
  line-height: 1.3rem;
  color: rgba(51, 51, 51, 0.5);
}
.form__button{
  width: 40%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
.form__addfile{
  width: 7rem;
  height: 7rem;
  position: relative;
  cursor: pointer;
}
.form__button--plusicon{
  width: 2rem;
  height: 2rem;
  object-fit: cover;
  position: absolute;
  top:  0;
  right: 0;
}
.form__button--fordericon{
  width: 7rem;
  height: 7rem;
  object-fit: cover;
}
.form__input{
  position: absolute;
  width: 0;
  height: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
}
.form__info{
  width: 70%;
  height: 60%;
  display: flex;
  margin: auto;
  background: #FFFFFF;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
}
.form__info--type{
  width: 30%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
.form__info--img{
  width: 5rem;
  height: 5rem;
  object-fit: cover;
}
.form__info--name,
.form__info--size{
  width: 60%;
  height: 100%;
  display: flex;
  align-items: center;
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.5rem;
  line-height: 1.5rem;
  color: rgb(51, 51, 51);
}
.form__info--blank{
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.6rem;
  line-height: 1.6rem;
  color: rgba(51, 51, 51, 0.5);
}
.form__button--apply{
  cursor: pointer;
  position: absolute;
  right: 3rem;
  width: 8rem;
  height: 80%;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  background: #3478F6;
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.5rem;
  line-height: 1.5rem;
  color: #FFFFFF;
}
</style>
