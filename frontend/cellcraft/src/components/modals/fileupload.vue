<template>
  <div>
    <form class="fileUpload-form" @submit.prevent="uploadFile">
        <label class="fileUpload-form__title">Select CSV to upload</label>
        <input df-file class="fileUpload-form__input" type="file" ref="selectFile" @change.prevent="previewFile" accept="text/csv" />
        <ul class="fileUpload-form__info" v-if="selectFile">
          <li>name : {{ selectFile.name }}</li>
          <li>size : {{ selectFile.size | formatBytes}}</li>
          <li>type : {{ selectFile.type }}</li>
        </ul>
        <input type="submit" value="업로드">
    </form>
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
      done_upload: false
    }
  },
  methods: {
    previewFile () {
      if (this.$refs.selectFile.files.length > 0) {
        this.selectFile = this.$refs.selectFile.files[0]
      }
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

<style>
.fileUpload-form{
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
}
.fileUpload-form__title{
    font-size: 1.5rem;
}

.fileUpload-form__input{
  width: 100px;
  margin: 2rem 0;
}

.fileUpload-form__info{
    margin-bottom: 1rem;
}
</style>
