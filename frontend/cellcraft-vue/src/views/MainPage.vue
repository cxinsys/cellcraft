<template>
  <div class="layout">
    <div class="main__bg-video">
      <video class="main-video" autoplay loop muted>
        <source src="@/assets/main_bg_fantastic.mp4" type="video/mp4" />
      </video>
    </div>
    <main class="main-component">
      <div class="main__textbox">
        <img class="main__name" :src="randomImage" alt="Random Image" />
        <span class="main__desc">
          Platform for Reconstruction of Inter-cell Type Gene Regulatory Network
          Model
          <br />
          based on Integrative Analysis of Single Cell RNA-seq, ATAC-seq and
          Spatial Transcriptomics Data
        </span>
        <div class="main__hotButtons">
          <button
            class="main__hotButtons-item"
            @click="scrollToBottom('menuId')"
          >
            Menu >
          </button>
          <button
            class="main__hotButtons-item"
            @click="scrollToBottom('introId')"
          >
            About >
          </button>
        </div>
      </div>
    </main>
    <div class="notice-component">
      <h1 class="notice__title">Notice</h1>
      <div class="notice__column">
        <p class="notice__info">
          [Notice] The lastest CellCraft release came out...
        </p>
        <p class="notice__date">2023.01.04</p>
      </div>
      <div class="notice__column">
        <p class="notice__info">[Event] Welcome to CellCraft workshop</p>
        <p class="notice__date">2023.01.04</p>
      </div>
      <img class="notice__list" src="@/assets/menu.png" />
    </div>
    <section class="menu-component" id="menuId">
      <ul class="menu-list">
        <li
          class="menu-list__card"
          @mouseover="mouseoverEvent"
          @mouseout="mouseoutEvent"
        >
          <img class="menu-list__icon" src="@/assets/e-book.png" />
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc"
              >The best tutorial you can see to get started with CellCraft</span
            >
          </div>
          <img class="arrow" src="@/assets/right-arrow.png" />
        </li>
        <li
          class="menu-list__card"
          @mouseover="mouseoverEvent"
          @mouseout="mouseoutEvent"
        >
          <img class="menu-list__icon" src="@/assets/document.png" />
          <div class="menu-list__textbox">
            <p class="menu-list__title">Document</p>
            <span class="menu-list__desc"
              >Algorithms, Data sets, Research results, detailed resources
              etc.</span
            >
          </div>
          <img class="arrow" src="@/assets/right-arrow.png" />
        </li>
        <li
          class="menu-list__card"
          @mouseover="mouseoverEvent"
          @mouseout="mouseoutEvent"
        >
          <img class="menu-list__icon" src="@/assets/chat.png" />
          <div class="menu-list__textbox">
            <p class="menu-list__title">Forum</p>
            <span class="menu-list__desc"
              >A community where users communicate freely with each other.</span
            >
          </div>
          <img class="arrow" src="@/assets/right-arrow.png" />
        </li>
        <li
          class="menu-list__card"
          @mouseover="mouseoverEvent"
          @mouseout="mouseoutEvent"
        >
          <img class="menu-list__icon" src="@/assets/database.png" />
          <div class="menu-list__textbox">
            <p class="menu-list__title">Dataset</p>
            <span class="menu-list__desc"
              >Single-cell data set accumulated through the study</span
            >
          </div>
          <img class="arrow" src="@/assets/right-arrow.png" />
        </li>
        <li
          class="menu-list__card"
          @mouseover="mouseoverEvent"
          @mouseout="mouseoutEvent"
        >
          <img class="menu-list__icon" src="@/assets/megaphone.png" />
          <div class="menu-list__textbox">
            <p class="menu-list__title">Notice</p>
            <span class="menu-list__desc"
              >Check out announcements, blogs, and events</span
            >
          </div>
          <img class="arrow" src="@/assets/right-arrow.png" />
        </li>
        <li
          class="menu-list__card"
          @mouseover="mouseoverEvent"
          @mouseout="mouseoutEvent"
        >
          <img class="menu-list__icon" src="@/assets/faq.png" />
          <div class="menu-list__textbox">
            <p class="menu-list__title">FAQ</p>
            <span class="menu-list__desc"
              >You can see frequently asked questions while using our
              service</span
            >
          </div>
          <img class="arrow" src="@/assets/right-arrow.png" />
        </li>
      </ul>
    </section>
    <div class="intro-component" id="introId">
      <div class="intro__textbox">
        <h1 class="intro__title">
          <span class="bold">ABOUT</span>
          <br />
          CELLCRAFT
        </h1>
        <ul class="intro__list">
          <li class="intro__list-item">Intro</li>
          <li class="intro__list-item">Pipeline</li>
          <li class="intro__list-item">License</li>
        </ul>
      </div>
      <div class="intro__descbox">
        <div class="intro__descbox__column">
          <h1 class="intro__desc-title">
            Single-cell analysis data platform for researchers
          </h1>
          <p class="intro__desc-text">
            CellCraft is a data platform for easy analysis of single cell data
            and efficient use of data. CellCraft contains many utilities for
            performing related tasks such as processing and quality control of
            high-throughput sequencing data and provides services on the web.
          </p>
          <p class="intro__desc-text">
            These tools are primarily designed for data analysts in
            bioinformatics, but even those who do not know programming can
            easily use CellCraft to handle data analysis. CellCraft was
            originally developed for single-cell analysis, but it is being
            developed to handle multiple biological data.
          </p>
        </div>
        <div class="intro__descbox__column">
          <img class="intro__desc-img" src="@/assets/example.jpg" />
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import logo1 from "@/assets/cellcraft_logo_text_5.png";
import logo2 from "@/assets/cellcraft_logo_text_8.png";

export default {
  methods: {
    mouseoverEvent(ev) {
      if (ev.target.children[2]) {
        ev.target.children[2].classList.add("right-move");
      }
    },
    mouseoutEvent(ev) {
      if (ev.target.children[2]) {
        ev.target.children[2].classList.remove("right-move");
      }
    },
    scrollToBottom(sectionId) {
      const section = document.querySelector(`#${sectionId}`);
      window.scrollTo({
        top: section.offsetTop,
        behavior: "smooth",
      });
    },
  },
  mounted() {
    const video = this.$refs.videoRef;
    const videoWidth = video.videoWidth;
    const videoHeight = video.videoHeight;
    const ratio = Math.max(
      window.innerWidth / videoWidth,
      window.innerHeight / videoHeight
    );
    video.style.width = videoWidth * ratio + "px";
    video.style.height = videoHeight * ratio + "px";
  },
  beforeDestroy() {
    const video = this.$refs.videoRef;
    video.pause();
  },
  data() {
    return {
      images: [logo1, logo2],
    };
  },
  computed: {
    randomImage() {
      const index = Math.floor(Math.random() * this.images.length);
      return this.images[index];
    },
  },
};
</script>

<style scoped>
.layout {
  width: 100%;
  /* height: 10%; */
  position: relative;
}
.bold {
  font-weight: 400;
  /* color: rgb(70, 70, 70); */
}
.arrow {
  display: inline-block;
  opacity: 0;
  width: 1.5rem;
  height: 1.5rem;
  object-fit: cover;
  position: absolute;
  top: 1rem;
  right: 3rem;
  transition: 0.5s;
}
.right-move {
  transform: translateX(2rem);
  opacity: 1;
}
/* .main-component, */
/* .notice-component, */
.menu-component {
  width: calc(100vw - 17px);
  max-width: calc(100% - 3rem);
  margin-right: auto;
  margin-left: auto;
  position: relative;
}
.main-component {
  top: 2rem;
  height: calc(100vh - 70px);
  background: linear-gradient(to bottom, #b7b7b75a, #000000d0);

  /* background: #000000ac; */
}
.notice-component {
  height: 70px;
  /* margin-top: 4rem; */
  padding: 0rem calc(50% - 43rem);
  display: flex;
  align-items: center;
  margin-right: auto;
  margin-left: auto;
  position: relative;
}
.notice-component:before {
  content: "";
  background: rgba(0, 85, 203, 0.8);
  width: calc(100vw - 17px);
  height: 100%;
  position: absolute;
  left: 0;
  /* margin-left: -100vw; */
  z-index: -1;
  top: 0;
}
.menu-component {
  padding: 6rem 0;
  display: flex;
  align-items: center;
  height: 55vh;
  box-sizing: border-box;
}
.intro-component {
  /* height: 22.5rem; */
  height: calc(45vh - 80px);
  display: flex;
  align-items: center;
  justify-content: center;
  margin-right: auto;
  margin-left: auto;
  position: relative;
}
.intro-component:before {
  content: "";
  background: linear-gradient(to bottom, #ffffffa8, #ffffff);
  width: calc(100vw - 17px);
  height: 100%;
  position: absolute;
  left: 0;
  /* margin-left: -100vw; */
  z-index: -1;
  top: 0;
}
.main__img {
  width: 90rem;
  max-width: calc(100%);
  height: inherit;
  object-fit: cover;
  object-position: left;
  border-radius: 1.5rem;
}
.main__bg-video {
  position: fixed;
  width: 100%;
  height: 100%;
  overflow: hidden;
  z-index: -1;
}
.main-video {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  object-fit: cover;
  /* width: 100rem; */
  /* max-width: calc(99%); */
  /* height: 20rem; */
  /* object-fit: cover;
  object-position: left; */
  /* border-radius: 1.5rem; */
}
.main__textbox {
  /* width: 60rem; */
  /* flex-grow: 1; */
  width: calc(100vw - 17px);
  max-width: calc(99%);
  height: 70rem;
  position: absolute;
  top: 0rem;
  left: 0;
  right: 0;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
}
.main__name {
  height: auto;
  /* object-fit: contain; */
  width: 60%;
  min-width: 30rem;
  max-width: 60rem;
}
.main__desc {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.04rem;
  line-height: 1.6rem;
  padding: 2rem 0;
  text-align: center;
  color: rgb(255, 255, 255);
}
.main__hotButtons {
  /* padding: 10rem; */
  display: flex;
  flex-direction: row;
  align-items: flex-end;
}
.main__hotButtons-item {
  border: none;
  background-color: transparent;
  color: rgb(109, 158, 235);
  font-size: 1.2rem;
  cursor: pointer;
  text-decoration: none;
  margin: 0.5rem 1rem;
}

.main__hotButtons-item:hover {
  text-decoration: underline;
}
.notice__title {
  width: 10%;
  height: 100%;
  display: flex;
  justify-content: flex-start;
  align-items: center;

  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 1.5rem;
  line-height: 1.5rem;
  text-align: left;
  color: #fff;
}
.notice__column {
  width: 45%;
  height: 100%;
  display: flex;
  align-items: center;
}
.notice__info,
.notice__date {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  text-align: center;
  color: #fff;
  opacity: 0.95;
  display: flex;
  align-items: center;
  justify-content: center;
}
.notice__info {
  width: 70%;
  height: 100%;
}
.notice__date {
  width: 30%;
  height: 100%;
}
.notice__list {
  width: 1.5rem;
  height: 1.5rem;
  object-fit: cover;
  filter: invert(91%) sepia(91%) saturate(30%) hue-rotate(198deg)
    brightness(106%) contrast(100%);
}
.menu-list {
  width: 100%;
  padding: 1rem 0;
  overflow: hidden;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  justify-content: center;
}
.menu-list__card {
  cursor: pointer;
  margin: 1rem 2rem;
  margin-bottom: 2rem;
  width: 24vw;
  max-width: calc(50% - 5rem);
  height: 16vh;
  display: flex;
  align-items: center;
  justify-content: center;
  background: #f3f3f3d2;

  /* box-shadow: rgba(0, 0, 0, 0.24) 0px 3px 8px; */
  box-shadow: rgba(0, 0, 0, 0.14) 0.5px 3px 5px;
  padding: 1rem;
  border-radius: 0.5rem;
  position: relative;
}
.menu-list__card:hover {
  box-shadow: rgba(0, 0, 0, 0.31) 0px 10px 20px, rgba(0, 0, 0, 0.35) 0px 6px 6px;
  transition: 0.5s;
  background: #f3f3f3;
}
.menu-list__icon {
  pointer-events: none;
  /* width: 2.5rem;
  height: 2.5rem; */
  width: 4rem;
  height: 4rem;
  object-fit: cover;
  margin-right: 2rem;
  opacity: 0.8;
}
.menu-list__textbox {
  pointer-events: none;
  width: 17.5rem;
  height: 5.5rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  position: relative;
}
.menu-list__title {
  pointer-events: none;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1.4rem;
  position: absolute;
  top: calc(100% - 5rem);
  /* color: rgb(70, 70, 70); */
  color: rgb(0, 85, 203);
}
.menu-list__desc {
  pointer-events: none;
  width: calc(100%);
  max-width: 19rem;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1.2rem;
  margin-left: 0.1rem;
  position: absolute;
  top: calc(100% - 3rem);
  color: rgb(80, 80, 80);
}
.intro__textbox {
  width: 20rem;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;
  border-right: 1px solid rgba(0, 0, 0, 0.3);
}
.intro__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 700;
  /* font-size: 2.75rem; */
  font-size: 2.3rem;
  line-height: 2.75rem;
  color: rgb(70, 70, 70);
  /* color: rgb(0, 85, 203); */
}
.intro__list {
  width: 100%;
  margin-top: 1rem;
  counter-reset: list-counter;
  list-style: none;
}
.intro__list-item {
  width: 100%;
  margin: 1rem 0;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 0.75rem;
  padding: 0.25rem;
}
.intro__list-item:hover {
  box-shadow: rgba(0, 0, 0, 0.4) -5px 5px, rgba(0, 0, 0, 0.3) -10px 10px;
  transform: translateX(0.5rem) translateY(-0.5rem);
  transition: 0.5s;
}
.intro__list-item:before {
  content: counter(list-counter);
  counter-increment: list-counter;
  width: 0.75rem;
  height: 0.75rem;
  padding: 0.5rem;
  margin-right: 0.5rem;
  border-radius: 50%;
  border: 0.25em solid #ccc;
  background: #000;
  color: #fff;
  font-family: arial;
  font-weight: bold;
  text-align: center;
  display: inline-block;
}
.intro__descbox {
  width: 70rem;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
.intro__descbox__column {
  width: 45%;
  height: 100%;
  margin: 0 1rem;
  /* margin-top: 10%; */
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}
.intro__descbox__column:last-child {
  /* margin-top: 1%; */
  align-items: center;
}
.intro__desc-title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 1.2rem;
  line-height: 1.4rem;
  color: rgb(50, 50, 50);
}
.intro__desc-text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 0.9rem;
  line-height: 1.5rem;
  margin-top: 1rem;
  opacity: 0.7;
  color: rgb(40, 40, 40);
}
.intro__desc-img {
  /* width: 28rem; */
  width: 28rem;
  height: 20rem;
  object-fit: cover;
  border-radius: 0.5rem;
}
</style>
