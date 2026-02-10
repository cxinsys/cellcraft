<template>
  <div class="layout" @click="ClickOut">
    <section class="folder">
      <div class="folder__header">
        <p class="folder__title left">Navigator</p>
        <div class="folder__createBtn right">
          <!-- <img class="folder__createBtn--icon" src="@/assets/add-folder.png" />
          <img class="folder__createBtn--icon" src="@/assets/refresh.png" /> -->
        </div>
      </div>
      <ul class="folder__list">
        <li class="folder__item" v-for="(folder, idx) in folders_list" :key="idx"
          v-bind:class="{ toggleFolder: toggleFolder === idx }" @click="folderClick(idx, folder[0])">
          <div class="folder__item--col">
            <img class="folder__item--icon" src="@/assets/arrow-bottom.png" v-if="toggleFolder === idx" />
            <img class="folder__item--icon" src="@/assets/arrow-right.png" v-else />
            <img class="folder__item--icon large" src="@/assets/open-folder.png" v-if="toggleFolder === idx" />
            <img class="folder__item--icon large" src="@/assets/folder.png" v-else />
          </div>
          <p class="folder__name">{{ folder[0] }}</p>
        </li>
      </ul>
    </section>
    <main class="files">
      <div class="files__header">
        <div class="header__column left">
          <p class="files__folder">{{ currentFolder }}</p>
        </div>
        <!-- <input
          class="files__search"
          type="text"
          name="search"
          placeholder="search anything..."
        /> -->
        <div class="header__column right">
          <label class="files__button">
            <img class="files__button--icon" src="@/assets/upload-file-black.png" />
            <h1>Upload File</h1>
            <input class="files__input" type="file" name="file" ref="selectFile" @change.prevent="uploadFile" />
          </label>
          <div class="progress__box" v-if="isUploading">
            <!-- Step 1: Upload progress bar -->
            <template v-if="uploadStep === 1">
              <div class="upload-bar-track">
                <div class="upload-bar-fill" :style="{ width: uploadPercentage + '%' }"></div>
              </div>
              <span class="upload-status-text">
                Uploading<template v-if="uploadTotalSteps > 1"> (1/{{ uploadTotalSteps }})</template>...
                {{ uploadPercentage }}%
              </span>
            </template>
            <!-- Step 2: Compression indeterminate bar (H5AD only) -->
            <template v-if="uploadStep === 2">
              <div class="upload-bar-track">
                <div class="upload-bar-fill upload-bar-fill--indeterminate"></div>
              </div>
              <span class="upload-status-text">Compressing (2/{{ uploadTotalSteps }})</span>
            </template>
          </div>
        </div>
      </div>

      <table class="files__table">
        <thead>
          <tr>
            <th>Name</th>
            <th>Date</th>
            <th>Type</th>
            <th>Size</th>
          </tr>
        </thead>
        <tbody>
          <tr class="files__item" v-for="(file, idx) in files_list" :key="idx" @contextmenu.prevent
            @click.right="RMouseClick($event, file.file_name, idx)" v-bind:class="{ select: R_Mouse_isActive }">
            <td>{{ extractFileName(file.file_name) }}</td>
            <td>{{ cutDateFromISO(file.created_at) }}</td>
            <td>{{ extractExtension(file.file_name) }}</td>
            <td>{{ formatBytes(file.file_size) }}</td>
          </tr>
        </tbody>
      </table>
      <ul class="toggle__menu" v-bind:class="{ open: R_Mouse_isActive }"
        :style="{ left: xPosition, top: yPosition }">
        <!-- <li>view</li>
        <li>plot</li>
        <li>rename</li> -->
        <li @click="removeFile">Delete</li>
      </ul>
    </main>
  </div>
</template>

<script>
import { uploadForm, getFiles, findFolder, deleteFile } from "@/api/index";
import { validateFileExtension } from "@/utils/validation";
import { generateUploadFileName } from "@/utils/filename";
import { calculateContextMenuPosition } from "@/utils/positionCalculator";
import { formatBytes, cutDateFromISO, extractFileName, extractExtension } from "@/utils/formatters";

export default {
  props: {
    doUploadFile: {
      type: String,
      default: "",
    },
  },

  data() {
    return {
      folders_list: [],
      files_list: [],
      R_Mouse_isActive: false,
      Clickout_isActive: false,
      toggleFolder: null,
      currentFolder: "data",
      xPosition: 0,
      yPosition: 0,
      selectFile: null,
      file_name: null,
      list_idx: null,
      uploadPercentage: 0,
      uploadStep: 0,        // 0=idle, 1=uploading, 2=compressing
      uploadTotalSteps: 1,  // 1=non-H5AD, 2=H5AD
      isDeletingFile: false, // Track file deletion state to prevent concurrent deletions
    };
  },

  computed: {
    isUploading() {
      return this.uploadStep > 0;
    },
  },

  methods: {
    // Formatter methods for template
    formatBytes,
    cutDateFromISO,
    extractFileName,
    extractExtension,

    RMouseClick(event, file_name, idx) {
      this.R_Mouse_isActive = false;
      const { x, y } = calculateContextMenuPosition(event.clientX, event.clientY);
      this.xPosition = x;
      this.yPosition = y;
      this.R_Mouse_isActive = true;
      this.file_name = file_name;
      this.list_idx = idx;
    },
    ClickOut() {
      this.R_Mouse_isActive = false;
    },
    folderClick(idx, folderName) {
      if (idx === this.toggleFolder) {
        this.toggleFolder = null;
      } else {
        this.toggleFolder = idx;
      }
      this.currentFolder = folderName;
    },
    async uploadFile() {
      if (this.$refs.selectFile.files.length > 0) {
        const file = this.$refs.selectFile.files[0];

        const validation = validateFileExtension(file.name);
        if (!validation.isValid) {
          alert(validation.message);
          return;
        }

        const newFileName = generateUploadFileName(this.currentFolder, file.name);
        this.selectFile = new File([file], newFileName);

        const isH5ad = file.name.toLowerCase().endsWith('.h5ad');

        // Initialize step state
        this.uploadTotalSteps = isH5ad ? 2 : 1;
        this.uploadStep = 1;
        this.uploadPercentage = 0;

        const form = new FormData();
        form.append("files", this.selectFile);

        const onUploadProgress = (progressEvent) => {
          if (!progressEvent.total || progressEvent.total === 0) return;
          const pct = Math.min(
            100,
            parseInt(Math.round((progressEvent.loaded * 100) / progressEvent.total))
          );
          if (pct >= 100 && isH5ad) {
            this.uploadPercentage = 100;
            this.$nextTick(() => {
              this.uploadStep = 2;
              this.uploadPercentage = 0;
            });
          } else {
            this.uploadPercentage = pct;
          }
        };

        try {
          await uploadForm(form, onUploadProgress);
          this.resetUploadState();

          const folderList = await findFolder({ folder_name: this.currentFolder });
          this.files_list = folderList.data;
        } catch (error) {
          console.error('File upload failed:', error);
          this.resetUploadState();

          let errorMessage = 'File upload failed. Please try again.';
          if (error.response) {
            if (error.response.data && error.response.data.detail) {
              errorMessage = error.response.data.detail;
            } else if (error.response.status === 413) {
              errorMessage = 'File is too large. Please choose a smaller file.';
            } else if (error.response.status === 400) {
              errorMessage = 'Invalid file. Please check the file format and try again.';
            } else if (error.response.status === 500) {
              errorMessage = 'Server error occurred. Please try again later.';
            }
          } else if (error.request) {
            errorMessage = 'Network error. Please check your connection.';
          }
          alert(errorMessage);
        }
      }
    },
    resetUploadState() {
      this.uploadStep = 0;
      this.uploadTotalSteps = 1;
      this.uploadPercentage = 0;
    },
    async removeFile() {
      // Prevent concurrent deletions
      if (this.isDeletingFile) {
        return;
      }

      // Show confirmation dialog BEFORE any state changes
      if (
        !confirm(
          "Are you sure you want to delete this file? This action cannot be undone."
        )
      ) {
        return; // User cancelled
      }

      // Set loading state
      this.isDeletingFile = true;

      try {
        const file = {
          file_name: this.file_name,
        };

        // Call API and wait for response
        await deleteFile(file);

        // Only remove from UI after successful deletion
        this.files_list.splice(this.list_idx, 1);

      } catch (error) {
        console.error('File deletion failed:', error);

        // User-friendly error messages
        let errorMessage = 'Failed to delete file. Please try again.';
        if (error.response?.data?.detail) {
          errorMessage = error.response.data.detail;
        } else if (error.response?.status === 404) {
          errorMessage = 'File not found. It may have been already deleted.';
        } else if (error.request) {
          errorMessage = 'Network error. Please check your connection.';
        }

        alert(errorMessage);
      } finally {
        // Reset loading state
        this.isDeletingFile = false;
      }
    },
    handleBeforeUnload(e) {
      if (this.uploadStep > 0) {
        e.preventDefault();
        e.returnValue = '';
      }
    },
  },

  beforeRouteLeave(to, from, next) {
    if (this.uploadStep > 0) {
      const answer = window.confirm(
        'File upload is in progress. Are you sure you want to leave this page?'
      );
      next(answer);
    } else {
      next();
    }
  },

  async mounted() {
    window.addEventListener('beforeunload', this.handleBeforeUnload);
    if (this.$route.query.doUploadFile) {
      this.$refs.selectFile.click();
    }
    try {
      const fileList = await getFiles();
      this.folders_list = fileList.data;
      const folderList = await findFolder({
        folder_name: this.currentFolder,
      });
      this.files_list = folderList.data;
    } catch (error) {
      console.error(error);
    }
  },

  beforeDestroy() {
    window.removeEventListener('beforeunload', this.handleBeforeUnload);
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  height: 100%;
  display: flex;
  position: relative;
  overflow: auto;
}

.folder {
  width: 18rem;
  height: 100%;
  border-right: 1px solid #aeaeae;
  background-color: rgb(201, 202, 203);
  /* width: 22rem;
  height: 100%;
  background: #cfcfcf;
  border-right: 1px solid #afafaf; */
}

.folder__header {
  width: 90%;
  height: 10%;
  margin: auto;
  display: flex;
  align-items: center;
  position: relative;
  color: rgba(0, 0, 0, 0.8);
}

.folder__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.7rem;
  line-height: 1.25rem;
  color: rgba(0, 0, 0, 0.8);
}

.folder__createBtn {
  width: 6rem;
  height: 2rem;
  vertical-align: middle;
}

.folder__createBtn--icon {
  width: 2rem;
  height: 2rem;
  object-fit: cover;
  margin: 0 0.5rem;
  opacity: 0.8;
}

.folder__list {
  width: 85%;
  height: 90%;
  margin: auto;
}

.folder__item {
  width: 100%;
  height: 5%;
  /* margin: 0 0rem; */
  border-radius: 0.5rem;
  display: flex;
  cursor: pointer;
}

.folder__item:hover {
  background: rgb(176, 177, 178);
}

.toggleFolder {
  background: rgb(176, 177, 178);
}

.folder__item--col {
  width: 5rem;
  height: 100%;
  display: flex;
  align-items: center;
}

.folder__item--icon {
  width: 1rem;
  height: 1rem;
  object-fit: contain;
  margin: 0 0.5rem;
}

.large {
  width: 1.5rem;
  height: 1.5rem;
}

.folder__name {
  display: flex;
  align-items: center;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
}

.files {
  width: 75%;
  height: 100%;
  background: #ffffff;
}

.files__header {
  width: 100%;
  height: 10%;
  display: flex;
  align-items: center;
  justify-content: center;
  position: relative;
}

.header__column {
  width: 50%;
  height: 100%;
  display: flex;
  align-items: center;
}

.right {
  position: absolute;
  right: 0;
}

.left {
  position: absolute;
  left: 0;
}

.files__folder {
  margin-left: 2rem;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.8rem;
  line-height: 1.25rem;
  color: rgba(0, 0, 0, 0.8);
  /* text-transform: capitalize; */
}

.files__button {
  width: 9rem;
  height: 2rem;
  padding: 0.2rem;
  display: flex;
  align-items: center;
  justify-content: center;
  border: none;
  /* background: #ffffff; */
  border-radius: 1.2rem;
  margin-right: 1rem;
  box-shadow: rgba(0, 0, 0, 0.15) 0px 0px 4px;
}

.files__button:hover {
  cursor: pointer;
  box-shadow: rgba(0, 0, 0, 0.35) 0px 0px 4px;
}

.files__button--icon {
  width: 1.75rem;
  height: 1.75rem;
  object-fit: contain;
  opacity: 0.8;
  margin-right: 0.5rem;
}

.files__input {
  display: none;
}

.files__search {
  width: 300px;
  height: 2.5rem;
  border: 1px solid #e1e1e1;
  border-radius: 1rem;
  padding: 0 2rem;
  outline-style: none;
  background: #f7f7f7;
}

.files__search:focus {
  border: 1px solid #bcbcbc;
}

.files__table {
  width: 95%;
  height: auto;
  margin: auto;
  border-collapse: collapse;
}

.files__table thead {
  height: 3rem;
  border-bottom: 1px solid #e1e1e1;
}

.files__table th {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.1rem;
  line-height: 1.1rem;

  text-align: left;
  vertical-align: middle;
  padding: 0 1rem;

  color: rgba(0, 0, 0, 0.5);
}

.files__table td {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  vertical-align: middle;
  padding: 1rem;
}

.files__item:hover {
  cursor: pointer;
  background: rgb(204, 218, 245);
}

.files__item td {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.files__header__img {
  margin-left: 5px;
  width: 15px;
}

.font-setting {
  font-size: 1vw;
  font-weight: bold;
  color: rgb(70, 70, 70);
}

.files__contents {
  width: 100vw;
  height: 75vh;

  margin-top: 20px;
}

.files__contents__item {
  display: flex;
  align-items: center;

  width: 100%;
  height: 40px;
}

/* .files__contents__item.select{
    background-color: #242F9B;
    color: white;
  } */

.files__contents__item:hover {
  background-color: #242f9b;
  color: white;
}

.toggle__menu {
  display: none;
  position: absolute;
  width: 200px;
  margin: 0;
  padding: 0;
  background: #ffffff;
  border-radius: 5px;
  list-style: none;
  box-shadow: 0 15px 35px rgba(50, 50, 90, 0.1), 0 5px 15px rgba(0, 0, 0, 0.07);
  overflow: hidden;
  z-index: 999999;
  text-transform: capitalize;
}

.toggle__menu.open {
  display: block;
  opacity: 1;
  position: absolute;
}

.toggle__menu>li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
}

.toggle__menu>li:hover {
  background: #e5e5e5;
}

.progress__box {
  display: flex;
  align-items: center;
  width: 18rem;
  height: 1.2rem;
  margin-right: 1rem;
}

.upload-bar-track {
  width: 10rem;
  height: 0.6rem;
  background-color: #eee;
  border-radius: 10px;
  overflow: hidden;
  flex-shrink: 0;
  margin-right: 0.5rem;
}

.upload-bar-fill {
  height: 100%;
  background-color: #4caf50;
  border-radius: 10px;
  transition: width 0.2s ease;
}

.upload-bar-fill--indeterminate {
  width: 40%;
  animation: indeterminate-slide 1.4s ease-in-out infinite;
}

@keyframes indeterminate-slide {
  0%   { transform: translateX(-100%); }
  50%  { transform: translateX(150%); }
  100% { transform: translateX(350%); }
}

.upload-status-text {
  font-family: "Montserrat", sans-serif;
  font-size: 0.8rem;
  color: #666;
  white-space: nowrap;
}
</style>
