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
            More
          </button>
          <!-- <button
            class="main__hotButtons-item"
            @click="scrollToBottom('introId')"
          >
            About >
          </button> -->
        </div>
      </div>
    </main>
    <!-- <div class="notice-component">
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
    </div> -->
    <div class="main-bump"></div>
    <section class="menu-component" id="menuId">
      <h1 class="menu__title">Tutorial1</h1>
      <h2 class="menu__description">
        Follow the posts below, ready the dataset for Cellcraft
      </h2>
      <ul class="menu-list">
        <div class="menu-list__col">
          <li
            class="menu-list__card"
            @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent"
          >
            <img class="menu-list__icon" src="@/assets/Seurat_work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">seurat</p>
              <span class="menu-list__desc"
                >A basic tutorial on single-cell transcriptomics using Seurat.
              </span>
            </div>
            <img class="arrow" src="@/assets/right-arrow.png" />
          </li>
          <li
            class="menu-list__card"
            @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent"
          >
            <img class="menu-list__icon" src="@/assets/scanpy_work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">scanpy</p>
              <span class="menu-list__desc"
                >A basic tutorial on single-cell transcriptomics using Scanpy.
              </span>
            </div>
            <img class="arrow" src="@/assets/right-arrow.png" />
          </li>
        </div>
        <div class="menu-list__col">
          <li
            class="menu-list__card"
            @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent"
          >
            <img class="menu-list__icon" src="@/assets/Monocle_Work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">monocle</p>
              <span class="menu-list__desc">
                Example of Trajectory Analysis Using Monocle
              </span>
            </div>
            <img class="arrow" src="@/assets/right-arrow.png" />
          </li>
          <li
            class="menu-list__card"
            @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent"
          >
            <img class="menu-list__icon" src="@/assets/PAGA_work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">PAGA</p>
              <span class="menu-list__desc"
                >Example of Trajectory Analysis Using PAGA</span
              >
            </div>
            <img class="arrow" src="@/assets/right-arrow.png" />
          </li>
        </div>
        <!-- <li
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
        </li> -->
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
import logo2 from "@/assets/cellcraft_logo_text_12.png";

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
  position: relative;
}

.bold {
  font-weight: 400;
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
  transform: translateX(1rem);
  opacity: 1;
}

.menu-component .main-component .main-bump .notice-component {
  width: 100%;
  max-width: 100%;
  margin-right: auto;
  margin-left: auto;
  position: relative;
}

.main-component {
  top: 2rem;
  height: calc(100vh - 0px);
  background-color: #000;
  animation-delay: 2s;
  animation: darkenAndLightenMain 10s ease forwards;
}

@keyframes darkenAndLightenMain {
  0% {
    background-color: #000;
  }
  100% {
    background-color: #000000ac;
  }
}

.main-bump {
  height: 30vh;
  animation: darkenAndLightenBump 10s ease forwards;
}

@keyframes darkenAndLightenBump {
  0% {
    background: linear-gradient(to bottom, #000000, #c3c3c3);
  }
  25% {
    background: linear-gradient(to bottom, #000000e8, #c3c3c3);
  }
  50% {
    background: linear-gradient(to bottom, #000000ca, #c3c3c3);
  }
  75% {
    background: linear-gradient(to bottom, #000000bc, #c3c3c3);
  }
  100% {
    background: linear-gradient(to bottom, #000000ac, #c3c3c3);
  }
}

.notice-component {
  height: 70px;
  padding: 0rem calc(50% - 43rem);
  display: flex;
  align-items: center;
  margin-right: auto;
  margin-left: auto;
  position: relative;
}

/* .notice-component:before {
  content: "";
  background: rgba(0, 85, 203, 0.8);
  width: calc(100vw - 0px);
  height: 100%;
  position: absolute;
  left: 0;
  z-index: -1;
  top: 0;
} */

.menu-component {
  padding: 6rem 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  height: 105vh;
  box-sizing: border-box;
  background-color: #c3c3c3;
}

.intro-component {
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
  background: linear-gradient(to bottom, #c3c3c3, #ffffff);
  width: 100%;
  height: 100%;
  position: absolute;
  left: 0;
  z-index: -1;
  top: 0;
}

.main__img {
  width: 90rem;
  max-width: 90%;
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
}

.main__textbox {
  width: 100%;
  max-width: 100%;
  height: calc(100vh - 70px);
  position: absolute;
  top: 35px;
  left: 0;
  right: 0;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
}

.main__name {
  height: auto;
  width: 40%;
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
  opacity: 0;
  animation-name: fadeIn;
  animation-duration: 6s;
  animation-fill-mode: forwards;
  animation-delay: 3s;
}

@keyframes fadeIn {
  from {
    opacity: 0;
  }
  to {
    opacity: 1;
  }
}

.main__hotButtons {
  display: flex;
  flex-direction: row;
  align-items: flex-end;
}

.main__hotButtons-item {
  border: none;
  background-color: transparent;
  color: rgb(143, 174, 225);
  font-size: 1.4rem;
  cursor: pointer;
  text-decoration: none;
  margin: 0.5rem 1rem;
  opacity: 0;
  animation-name: fadeIn;
  animation-duration: 3s;
  animation-fill-mode: forwards;
  animation-delay: 6s;
}

.main__hotButtons-item:hover {
  text-decoration: underline;
}

/* .notice__title {
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
} */

/* .notice__info,
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
} */

.menu__title {
  font-family: "Montserrat", sans-serif; /* Montserrat 폰트 적용 */
  font-weight: 700; /* 두께는 bold로 설정 */
  font-size: 2.5rem; /* 크기는 2.5rem으로 설정 */
  color: #333333; /* 글자 색상 설정 */
  text-transform: uppercase; /* 모든 글자를 대문자로 */
  letter-spacing: 0.05em; /* 글자 사이 간격 설정 */
  text-align: center; /* 중앙 정렬 */
  margin-top: 10px; /* 상단 여백 */
  margin-bottom: 10px; /* 하단 여백 */
}

.menu__description {
  font-family: "Montserrat", sans-serif; /* Montserrat 폰트 적용 */
  font-weight: 700; /* 두께는 bold로 설정 */
  font-size: 1.9rem; /* 크기는 2.5rem으로 설정 */
  color: #525252; /* 글자 색상 설정 */
  text-transform: uppercase; /* 모든 글자를 대문자로 */
  letter-spacing: 0.05em; /* 글자 사이 간격 설정 */
  text-align: center; /* 중앙 정렬 */
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

.menu-list__col {
  width: 40%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: space-evenly;
  position: relative;
}

/* .menu-list__col:first-child::before,
.menu-list__col:first-child::after {
  content: "▼";
  position: absolute;
  bottom: -25px;
  font-size: 2rem;
  color: #ffffff;
  animation: moveArrow 1.5s infinite;
} */
/* 
.menu-list__col:first-child::before {
  left: calc(50% - 23.5rem);
}

.menu-list__col:first-child::after {
  left: calc(50% + 23.5rem);
} */

.menu-list__card {
  cursor: pointer;
  margin: 2rem 1rem;
  margin-bottom: 2rem;
  width: 34rem;
  height: 13rem;
  max-width: 80%;
  max-height: 80%;
  display: flex;
  align-items: center;
  justify-content: center;
  background: #f3f3f3d2;
  box-shadow: rgba(0, 0, 0, 0.34) 0px 0px 10px;
  padding: 1rem;
  border-radius: 2.5rem;
  /* transition: 0.5s; */
  position: relative;
}

.menu-list__card:hover {
  width: 36rem;
  max-width: calc(80% + 2rem);
  box-shadow: rgba(0, 0, 0, 0.31) 0px 5px 20px, rgba(0, 0, 0, 0.35) 0px 3px 10px;
  transition: 0.5s;
  background: #f3f3f3;
}

.menu-list__card:first-child::before {
  content: "▼";
  position: absolute;
  left: calc(50% - 1rem);
  font-size: 2rem;
  color: #ffffff;
  animation: moveArrow 3s infinite;
}

.menu-list__icon {
  pointer-events: none;
  width: 100%;
  height: 100%;
  object-fit: cover;
  position: absolute;
  opacity: 0.25;
  border-radius: 2.5rem;
}

.menu-list__textbox {
  pointer-events: none;
  width: 17.5rem;
  height: 5.5rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  position: relative;
  border-radius: 0.5rem;
}

.menu-list__title {
  pointer-events: none;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1.4rem;
  position: absolute;
  top: 0rem;
  color: rgb(31, 56, 90);
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
  margin-left: 2rem;
}

.intro__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 700;
  font-size: 2.3rem;
  line-height: 2.75rem;
  color: rgb(70, 70, 70);
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
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.intro__descbox__column:last-child {
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
  width: 28rem;
  height: 20rem;
  object-fit: cover;
  border-radius: 0.5rem;
}

/* Defining the animation */
@keyframes moveArrow {
  0%,
  100% {
    bottom: -51px; /* Starting and ending position */
  }
  50% {
    bottom: -44px; /* Middle position */
  }
}
</style>
