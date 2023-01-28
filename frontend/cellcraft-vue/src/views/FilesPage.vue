<style scoped>
.files__header {
  width: 100vw;
  height: 5vh;
  display: flex;
  align-items: center;
  margin-top: 50px;
  border-bottom: 1px #a1a1a1 solid;
}
.l-col {
  width: 25%;
  height: 100%;

  display: flex;
  align-items: center;
  justify-content: center;
}
.l-col:first-child {
  width: 50%;
  height: 100%;

  display: flex;
  align-items: center;
  justify-content: center;
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

.files_menu {
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
}
.files_menu.open {
  display: block;
  opacity: 1;
  position: absolute;
}
.files_menu > li {
  border-left: 3px solid transparent;
  transition: ease 0.2s;
  padding: 10px;
}
.files_menu > li:hover {
  background: #e5e5e5;
}
</style>

<template>
  <section>
    <div class="files__header">
      <div class="l-col">
        <div class="font-setting">Name</div>
        <img
          class="files__header__img"
          src="@/assets/descending-filter-desc.png"
          alt=""
        />
      </div>
      <div class="l-col">
        <div class="font-setting">Size</div>
      </div>
      <div class="l-col">
        <div class="font-setting">Date</div>
      </div>
    </div>

    <ul class="files__contents" @click="ClickOut">
      <li
        class="files__contents__item"
        v-for="(file, idx) in files_list"
        :key="idx"
        @contextmenu.prevent
        @click.right="RMouseClick($event)"
        v-bind:class="{ select: R_Mouse_isActive }"
      >
        <div class="l-col">{{ file.file_name }}</div>
        <div class="l-col">{{ file.file_size }}</div>
        <div class="l-col">{{ file.created_at }}</div>
      </li>
      <ul
        ref="filesMenu"
        class="files_menu"
        v-bind:class="{ open: R_Mouse_isActive }"
        :style="{ left: xPosition, top: yPosition }"
      >
        <li>view</li>
        <li>plot</li>
        <li>rename</li>
        <li>delete</li>
      </ul>
    </ul>
  </section>
</template>

<script>
import { getFiles } from "@/api/index";

export default {
  data() {
    return {
      files_list: [],
      R_Mouse_isActive: false,
      Clickout_isActive: false,
    };
  },

  methods: {
    // mousePos (event) {

    //   console.log(this.xPosition)
    // },
    async RMouseClick(event) {
      this.R_Mouse_isActive = false;
      this.xPosition = event.clientX + "px";
      this.yPosition = event.clientY + "px";
      this.R_Mouse_isActive = true;
    },
    ClickOut() {
      this.R_Mouse_isActive = false;
    },
  },

  async mounted() {
    const fileList = await getFiles();
    console.log(fileList.data);
    this.files_list = fileList.data;
  },
};
</script>
