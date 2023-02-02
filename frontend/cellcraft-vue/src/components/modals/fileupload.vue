<template>
  <div class="layout">
    <form class="cloud-form" @submit.prevent="uploadFile">
      <div class="cloud-row">
        <label class="form__button">
          <div class="form__addfile">
            <img class="form__button--foldericon" src="@/assets/add-file.png" />
            <!-- <img class="form__button--plusicon" src="@/assets/plus.png" /> -->
          </div>
          <input
            df-file
            class="form__input"
            type="file"
            ref="selectFile"
            @change.prevent="previewFile"
            accept="text/csv"
          />
        </label>
        <div class="cloud-textbox">
          <h1 class="cloud-textbox__title">Get file</h1>
          <p class="cloud-textbox__desc">
            Get file from CELLCRAFT files directory.
          </p>
        </div>
      </div>
      <div class="cloud-row">
        <label class="form__button">
          <div class="form__addfile">
            <img
              class="form__button--foldericon"
              src="@/assets/upload-file.png"
            />
            <!-- <img class="form__button--plusicon" src="@/assets/plus.png" /> -->
          </div>
          <input
            df-file
            class="form__input"
            type="file"
            ref="selectFile"
            @change.prevent="previewFile"
            accept="text/csv"
          />
        </label>
        <div class="cloud-textbox">
          <h1 class="cloud-textbox__title">Upload file</h1>
          <p class="cloud-textbox__desc">
            Upload file to CELLCRAFT files directory.
          </p>
        </div>
      </div>
      <div class="cloud-row">
        <label class="form__button">
          <div class="form__addfile">
            <img
              class="form__button--foldericon"
              src="@/assets/manage-file.png"
            />
            <!-- <img class="form__button--plusicon" src="@/assets/plus.png" /> -->
          </div>
          <input
            df-file
            class="form__input"
            type="file"
            ref="selectFile"
            @change.prevent="previewFile"
            accept="text/csv"
          />
        </label>
        <div class="cloud-textbox">
          <h1 class="cloud-textbox__title">Files directory</h1>
          <p class="cloud-textbox__desc">Open CELLCRAFT files directory.</p>
        </div>
      </div>
    </form>
    <form class="fileUpload-form">
      <div class="form-row">
        <h1 class="form__name">Choose recent file</h1>
        <div class="form__selectFile">
          <ul class="form__fileList">
            <li
              class="fileList__item"
              v-for="(file, idx) in files_list"
              :key="idx"
              @click="fileSelect"
            >
              <p class="fileList__text">{{ file.file_name }}</p>
            </li>
            <li class="fileList__item" v-if="files_list.length === 0">
              <p class="fileList__text--blank">Please upload a new file</p>
            </li>
          </ul>
        </div>
      </div>
      <div class="form-row">
        <h1 class="form__name">Current file</h1>
      </div>
      <div class="form-row">
        <ul class="form__info" v-if="selectFile">
          <li class="form__info--type">
            <img class="form__info--img" src="@/assets/csv.png" />
          </li>
          <li class="form__info--name">
            {{ selectFile.name }}&nbsp;&nbsp;&nbsp;{{
              selectFile.size | formatBytes
            }}
          </li>
        </ul>
        <ul class="form__info" v-else>
          <li class="form__info--blank">Please add data file</li>
        </ul>
        <label class="form__button--apply">
          Apply
          <input class="form__input" type="submit" value="업로드" />
        </label>
      </div>
    </form>
  </div>
</template>

<script>
import { uploadForm, getFiles, findFile } from "@/api/index";

export default {
  data() {
    return {
      node_name: "File",
      selectFile: null,
      is_upload: false,
      done_upload: false,
      files_list: [],
    };
  },
  methods: {
    previewFile() {
      if (this.$refs.selectFile.files.length > 0) {
        this.selectFile = this.$refs.selectFile.files[0];
      }
    },
    async uploadFile() {
      if (this.selectFile === this.$refs.selectFile.files[0]) {
        this.is_upload = true;
        const form = new FormData();
        form.append("files", this.selectFile);
        const response = await uploadForm(form);
        console.log(response);
        if (response) {
          this.done_upload = true;
          const fileList = await getFiles();
          console.log(fileList.data);
          this.files_list = fileList.data;
        }
      }
    },
    async fileSelect(ev) {
      console.dir(ev.target.innerText);
      const fileInfo = {
        file_name: ev.target.innerText,
      };
      const selectFile = await findFile(fileInfo);
      this.selectFile = {
        name: selectFile.data.file_name,
        size: selectFile.data.file_size,
      };
      console.log(selectFile.data);
    },
  },
  filters: {
    formatBytes(a, b) {
      if (a === 0) return "0 Bytes";
      const c = 1024;
      const d = b || 2;
      const e = ["Bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"];
      const f = Math.floor(Math.log(a) / Math.log(c));

      return parseFloat((a / Math.pow(c, f)).toFixed(d)) + " " + e[f];
    },
  },
  async mounted() {
    const fileList = await getFiles();
    console.log(fileList.data);
    this.files_list = fileList.data;
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  position: relative;
  display: flex;
  align-items: flex-start;
}
.cloud-form {
  width: 45%;
  height: 95%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  margin: 1rem 0 1rem 1rem;
  padding: 3rem 0 3rem 1rem;
  border-radius: 1rem;
  box-sizing: border-box;
  background-color: rgb(255, 255, 255);
}
.cloud-row {
  width: 100%;
  padding: 5% 0;
  display: flex;
  flex-direction: row;
  align-items: flex-start;
}
.cloud-textbox {
  width: 17rem;
  height: 3rem;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;
  padding: 0 0 0 1rem;
}
.cloud-textbox__title {
  color: rgb(51, 51, 51);
  font-size: 1rem;
  font-weight: 400;
}
.cloud-textbox__desc {
  color: rgb(51, 51, 51);
  font-size: 1rem;
  font-weight: 200;
  padding: 0.2rem 0 0 0;
}
.fileUpload-form {
  width: 55%;
  height: 95%;
  margin: 1rem 2rem 1rem 0;
  display: flex;
  flex-direction: column;
  align-items: center;
}
.form-row {
  width: 100%;
  height: 8%;
  position: relative;
}
.form-row:first-child {
  height: 75%;
  margin: 1rem;
}
.form-row:last-child {
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  height: 17%;
  padding: 0 2rem;
  box-sizing: border-box;
}
.form__name {
  width: 100%;
  height: 3rem;
  display: flex;
  align-items: center;
  padding: 0 2rem;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.4rem;
  line-height: 1.4rem;
  color: rgb(51, 51, 51);
}
.form__selectFile {
  width: 100%;
  height: 85%;
  display: flex;
  padding: 0 0 0 0.5rem;
}
.form__fileList {
  width: 80%;
  height: 100%;
  display: flex;
  padding: 0% 3%;
  flex-direction: column;
  align-items: center;
  overflow: hidden;
}
.fileList__item {
  width: 90%;
  padding: 0.5rem;
  background: #ffffff;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  margin-bottom: 0.1rem;
}
.fileList__text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgb(53, 54, 58);
}
.fileList__text--blank {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgba(51, 51, 51, 0.5);
}
.form__button {
  width: 3rem;
  height: 3rem;
  display: flex;
  align-items: center;
  justify-content: center;
}
.form__addfile {
  width: 3rem;
  height: 3rem;
  position: relative;
  cursor: pointer;
  background-color: rgb(241, 243, 244);
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 1rem;
}
.form__button--plusicon {
  width: 1rem;
  height: 1rem;
  object-fit: cover;
  position: absolute;
  top: 2rem;
  right: 2rem;
}
.form__button--foldericon {
  width: 70%;
  height: 70%;
  object-fit: contain;
  position: absolute;
  top: 0.5rem;
  right: 0.5rem;
}
.form__input {
  position: absolute;
  width: 0;
  height: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
}
.form__info {
  width: 80%;
  height: 60%;
  display: flex;
  /* margin: auto; */
  top: 0;
  background: #ffffff;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
}
.form__info--type {
  width: 30%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
.form__info--img {
  width: 2rem;
  height: 2rem;
  object-fit: cover;
}
.form__info--name,
.form__info--size {
  width: 60%;
  height: 100%;
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgb(51, 51, 51);
}
.form__info--blank {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  color: rgba(51, 51, 51, 0.5);
}
.form__button--apply {
  cursor: pointer;
  position: absolute;
  left: 80%;
  top: 7%;
  width: 20%;
  height: 46%;
  box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
  border-radius: 0.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
  background: rgb(170, 193, 240);
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  color: rgb(244, 246, 251);
}

@media (prefers-color-scheme: dark) {
  .form__button--foldericon {
    filter: invert(97%) sepia(99%) saturate(0%) hue-rotate(123deg)
      brightness(107%) contrast(101%);
  }
  .cloud-form {
    width: 45%;
    height: 95%;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    margin: 1rem 0 1rem 1rem;
    padding: 3rem 0 3rem 1rem;
    border-radius: 1rem;
    box-sizing: border-box;
    background-color: rgb(41, 43, 48);
  }
  .form__addfile {
    width: 3rem;
    height: 3rem;
    position: relative;
    cursor: pointer;
    background-color: rgb(32, 34, 39);
    box-shadow: 0px 4px 4px rgba(0, 0, 0, 0.25);
    border-radius: 1rem;
  }
  .cloud-textbox__title {
    color: rgb(255, 255, 255);
    font-size: 1rem;
    font-weight: 400;
  }
  .cloud-textbox__desc {
    color: rgb(255, 255, 255);
    font-size: 1rem;
    font-weight: 200;
    padding: 0.2rem 0 0 0;
  }
  .form__name {
    width: 100%;
    height: 3rem;
    display: flex;
    align-items: center;
    padding: 0 2rem;
    font-family: "Montserrat", sans-serif;
    font-style: normal;
    font-weight: 400;
    font-size: 1.4rem;
    line-height: 1.4rem;
    color: rgb(255, 255, 255);
  }
}
</style>
