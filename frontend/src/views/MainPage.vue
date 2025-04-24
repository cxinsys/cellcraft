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
          Application for reconstructing gene regulatory networks (GRNs)
          <br />
          based on single-cell RNA sequencing (scRNA-seq) data
        </span>


        <div class="main__hotButtons">
          <button class="main__hotButtons-item" @click="scrollToBottom('introId')">
            More
          </button>
        </div>
      </div>
    </main>
    <div class="main-bump"></div>
    <div class="intro-component" id="introId">
      <div class="intro__descbox">
        <div class="intro__descbox__column">
          <h1 class="menu__title">{{ textDict.intro.title }}</h1>
          <p class="intro__desc-text">
            {{ textDict.intro.desc }}
          </p>
          <p v-for="(line, index) in textDict.intro.desc1" :key="'intro-desc-' + index" class="intro__desc-text">
            {{ line }}
          </p>
        </div>
        <div class="intro__descbox__column">
          <h1 class="menu__title">{{ textDict.gettingStarted.title }}</h1>
          <p class="intro__desc-text">
            {{ textDict.gettingStarted.desc }}
          </p>
          <p v-for="(line, index) in textDict.gettingStarted.desc1" :key="'getting-desc-' + index"
            class="intro__desc-text">
            {{ line }}
          </p>
        </div>
      </div>
    </div>
    <section class="menu-component" id="menuId">
      <h1 class="menu__title">Tutorial</h1>
      <h2 class="menu__description">
        Follow the posts below, ready the dataset for CellCraft
      </h2>
      <ul class="menu-list">
        <div class="menu-list__col">
          <li ref="seurat" class="menu-list__card" @click="openHtmlPage('seurat')" @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent">
            <img class="menu-list__icon" src="@/assets/Seurat_work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">Seurat</p>
              <span class="menu-list__desc">A basic tutorial on single-cell transcriptomics using Seurat.
              </span>
            </div>
          </li>
          <li ref="monocle" class="menu-list__card" @click="openHtmlPage('monocle')" @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent">
            <img class="menu-list__icon" src="@/assets/Monocle_Work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">Monocle</p>
              <span class="menu-list__desc">
                A tutorial on trajectory analysis using Monocle.
              </span>
            </div>
          </li>
        </div>
        <div class="menu-list__col">
          <li ref="scanpy" @click="openHtmlPage('scanpy')" class="menu-list__card" @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent">
            <img class="menu-list__icon" src="@/assets/scanpy_work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">Scanpy</p>
              <span class="menu-list__desc">A basic tutorial on single-cell transcriptomics using Scanpy.
              </span>
            </div>
          </li>
          <li ref="PAGA" @click="openHtmlPage('PAGA')" class="menu-list__card" @mouseover="mouseoverEvent"
            @mouseout="mouseoutEvent">
            <img class="menu-list__icon" src="@/assets/PAGA_work.png" />
            <div class="menu-list__textbox">
              <p class="menu-list__title">PAGA</p>
              <span class="menu-list__desc">A tutorial on trajectory analysis using PAGA.</span>
            </div>
          </li>
        </div>
      </ul>
    </section>
    <section class="menu-component">
      <div class="sub-description__box">
        <h2 class="menu__description sub-description">Convert Output</h2>
        <h2 class="menu__description sub-description">multi omics</h2>
      </div>
      <ul class="menu-list">
        <div class="menu-list__col">
          <li ref="output-conversion" @click="openHtmlPage('output-conversion')" class="menu-list__largeCard"
            @mouseover="mouseoverEvent" @mouseout="mouseoutEvent">
            <img class="menu-list__largeIcon" src="@/assets/output-conversion.png" />
            <div class="menu-list__textbox largeTextbox">
              <p class="menu-list__title">
                {{ textDict.outputConversion.title }}
              </p>
              <span class="menu-list__desc">{{ textDict.outputConversion.desc }}
              </span>
            </div>
          </li>
        </div>
        <div class="menu-list__col">
          <li ref="multi-omics" @click="openHtmlPage('multi-omics')" class="menu-list__largeCard"
            @mouseover="mouseoverEvent" @mouseout="mouseoutEvent">
            <img class="menu-list__largeIcon" src="@/assets/multi-omics.png" />
            <div class="menu-list__textbox largeTextbox">
              <p class="menu-list__title">
                {{ textDict.multiOmics.title }}
              </p>
              <span class="menu-list__desc">{{ textDict.multiOmics.desc }}
              </span>
            </div>
          </li>
        </div>
      </ul>
    </section>
    <section class="about-section">
      <div class="about-title-block">
        <h2 class="about-subtitle">ABOUT</h2>
        <h1 class="about-title">CELLCRAFT</h1>
        <p class="about-tagline">Single-cell GRN analysis without coding</p>
      </div>

      <div class="about-container">
        <div class="about-content-box">
          <div class="about-description">
            <p><strong>CellCraft</strong> is a web-based application for reconstructing gene regulatory networks (GRNs)
              from single-cell RNA sequencing (scRNA-seq) data.</p>
            <p>It features an intuitive visual interface that integrates seven GRN inference tools—including
              <strong>TENET</strong> and <strong>FastTENET</strong>, developed by our research team—through modular
              plugin integration.
            </p>
            <p>CellCraft enables researchers to <strong>configure workflows</strong>, <strong>run analysis</strong>, and
              <strong>explore results visually</strong>, all without needing programming knowledge.
            </p>
          </div>
        </div>
      </div>

      <div class="about-links-centered">
        <a href="https://github.com/cxinsys/cellcraft" target="_blank" class="link-card">
          <img src="@/assets/github.png" alt="GitHub" />
          <span>GitHub</span>
        </a>
        <a href="https://cellcraft.gitbook.io/cellcraft-docs" target="_blank" class="link-card">
          <img src="@/assets/gitbook.png" alt="Docs" />
          <span>Docs</span>
        </a>
        <a href="https://www.youtube.com/@CellCraft-cislab" target="_blank" class="link-card">
          <img src="@/assets/youtube.png" alt="YouTube" />
          <span>YouTube</span>
        </a>
      </div>
    </section>

  </div>
</template>

<script>
import logo1 from "@/assets/cellcraft_logo_text_5.png";
import logo2 from "@/assets/cellcraft_logo_text_12.png";

export default {
  data() {
    return {
      images: [logo1, logo2],
      textDict: {
        "seurat": {
          "title": "Seurat",
          "desc": "A basic tutorial on single-cell transcriptomics using Seurat."
        },
        "monocle": {
          "title": "Monocle",
          "desc": "A tutorial on trajectory analysis using Monocle."
        },
        "scanpy": {
          "title": "Scanpy",
          "desc": "A basic tutorial on single-cell transcriptomics using Scanpy."
        },
        "PAGA": {
          "title": "PAGA",
          "desc": "A tutorial on trajectory analysis using PAGA."
        },
        "outputConversion": {
          "title": "Output Conversion",
          "desc": "A guide to converting analyzed files (RDS to H5AD)."
        },
        "multiOmics": {
          "title": "Multi-Omics",
          "desc": "A tutorial on multi-omics analysis (scRNA-seq + scATAC-seq) using Seurat."
        },
        intro: {
          title: "Intro",
          desc: "Welcome to the official tutorial for CellCraft.",
          desc1: [
            "CellCraft offers a user-friendly interface that allows researchers to quickly explore, analyze, and visualize gene regulatory networks.",
            "CellCraft can process complex datasets, providing accurate and actionable insights into genetic interactions.",
            "Users can upload their data, select parameters for analysis, and witness the transformative power of CellCraft and TENET.",
            "This main tutorial aims to give users a comprehensive overview of CellCraft, guiding them through its features and capabilities.",
          ],
        },
        gettingStarted: {
          title: "Getting Started",
          desc: "Navigating CellCraft is simple and intuitive.",
          desc1: [
            "Refer to the sub-pages for detailed guidance on each step of the workflow—from dataset upload to GRN inference and result visualization.",
            "Whether conducting a quick test or a more complex multi-step analysis, CellCraft helps researchers efficiently carry out GRN studies using integrated tools and visual programming."
          ],
        }
      },
    };
  },
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
    openHtmlPage(ref) {
      const htmlName = ref;
      const routePath = this.$router.resolve({
        path: "/tutorial",
        query: { htmlName },
      }).href;
      window.open(routePath, "_blank");
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
  computed: {
    randomImage() {
      const index = Math.floor(Math.random() * this.images.length);
      return this.images[index];
    },
  },
};
</script>

<style scoped>
:root {
  --primary-color: #1f385a;
  --secondary-color: #f4f7fc;
  --text-color: #333;
  --text-light: #666;
  --text-white: #fff;
  --background-light: #f3f3f3;
  --shadow-color: rgba(0, 0, 0, 0.1);
  --transition-speed: 0.3s;
}

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

.main-bump {
  height: 30vh;
  animation: darkenAndLightenBump 10s ease forwards;
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

.menu-component {
  padding: 5rem 0 0rem 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  /* height: 105vh; */
  box-sizing: border-box;
  background-color: #c3c3c3;
}

.intro-component,
.about-component {
  /* height: calc(35vh - 80px); */
  display: flex;
  align-items: center;
  justify-content: center;
  margin-right: auto;
  margin-left: auto;
  position: relative;
  padding-top: 5rem;
}

.intro-component {
  background-color: #c3c3c3;
}

.about-component::before {
  content: "";
  background: linear-gradient(to bottom, #c3c3c3, #ffffff);
  width: 100%;
  height: 100%;
  position: absolute;
  left: 0;
  z-index: -1;
  top: 0;
}

/* .intro__descbox {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
} */

/* .intro__descbox__column {
  width: 70%;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;
  padding: 9.7rem 0 2rem 0;
} */

.intro__descbox__column {
  border-bottom: 1px solid rgba(0, 0, 0, 0.3);
}

.intro__desc-text {
  margin-left: 0.2rem;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1.1rem;
  line-height: 1.5rem;
  text-align: center;
  color: rgb(82, 82, 82);
}

.intro__descbox {
  display: flex;
  flex-wrap: wrap;
  justify-content: center;
  align-items: stretch;
  gap: 2rem;
  padding: 5rem 2rem;
}

.intro__descbox__column {
  flex: 1 1 500px;
  max-width: 600px;
  min-height: 100%;
  /* 높이 강제 통일을 위한 설정 */
  padding: 2rem;
  display: flex;
  flex-direction: column;
  justify-content: flex-start;
  background-color: #fff;
  border-radius: 1rem;
  box-shadow: 0 4px 16px rgba(0, 0, 0, 0.1);
}

.menu__title {
  font-size: 1.5rem;
  font-weight: 700;
  color: #222;
  margin-bottom: 1rem;
  text-align: left;
}

.intro__desc-text {
  text-align: left;
  color: #444;
  font-size: 1rem;
  line-height: 1.5rem;
  margin-bottom: 1rem;
  /* 문장 사이 간격 */
}



.main__img {
  width: 90rem;
  max-width: 90%;
  height: inherit;
  object-fit: cover;
  object-position: left;
  border-radius: 1rem;
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
  width: 100%;
  font-family: "Montserrat", sans-serif;
  /* Montserrat 폰트 적용 */
  font-weight: 700;
  /* 두께는 bold로 설정 */
  font-size: 2.5rem;
  /* 크기는 2.5rem으로 설정 */
  color: #333333;
  /* 글자 색상 설정 */
  text-transform: uppercase;
  /* 모든 글자를 대문자로 */
  letter-spacing: 0.05em;
  /* 글자 사이 간격 설정 */
  text-align: center;
  /* 중앙 정렬 */
  margin-top: 10px;
  /* 상단 여백 */
  margin-bottom: 1rem;
  /* 하단 여백 */
}

.menu__description {
  font-family: "Montserrat", sans-serif;
  /* Montserrat 폰트 적용 */
  font-weight: 700;
  /* 두께는 bold로 설정 */
  font-size: 1.9rem;
  /* 크기는 2.5rem으로 설정 */
  color: #525252;
  /* 글자 색상 설정 */
  text-transform: uppercase;
  /* 모든 글자를 대문자로 */
  letter-spacing: 0.05em;
  /* 글자 사이 간격 설정 */
  text-align: center;
  /* 중앙 정렬 */
  width: 70%;
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
  width: 26rem;
  height: 15rem;
  max-width: 80%;
  max-height: 80%;
  display: flex;
  align-items: center;
  justify-content: center;
  background: #f3f3f3d2;
  box-shadow: rgba(0, 0, 0, 0.34) 0px 0px 10px;
  padding: 1rem;
  border-radius: 1rem;
  /* transition: 0.5s; */
  position: relative;
}

.menu-list__card:hover {
  /* width: 36rem;
  max-width: calc(80% + 2rem); */
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
  bottom: -50px;
  /* Middle position */
}

.menu-list__icon {
  /* pointer-events: none; */
  width: 100%;
  height: 100%;
  object-fit: cover;
  position: absolute;
  opacity: 0.05;
  border-radius: 1rem;
  transition: opacity 1.5s ease;
  /* 투명도가 서서히 변화하는 트랜지션 설정 */
}

.menu-list__icon:hover {
  opacity: 0.5;
}

.menu-list__textbox {
  pointer-events: none;
  width: 17.5rem;
  height: 5.5rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  text-align: center;
  position: relative;
  border-radius: 1rem;
}

.menu-list__title {
  pointer-events: none;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  line-height: 1.4rem;
  text-align: center;
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

.about__textbox {
  width: 20rem;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;
  border-right: 1px solid rgba(0, 0, 0, 0.3);
  margin-left: 2rem;
}

.about__title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 700;
  font-size: 2.3rem;
  line-height: 2.75rem;
  color: rgb(70, 70, 70);
}

.about__list {
  width: 100%;
  margin-top: 1rem;
  counter-reset: list-counter;
  list-style: none;
}

.about__list-item {
  width: 100%;
  margin: 1rem 0;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 0.75rem;
  padding: 0.25rem;
}

.about__list-item:hover {
  box-shadow: rgba(0, 0, 0, 0.4) -5px 5px, rgba(0, 0, 0, 0.3) -10px 10px;
  transform: translateX(0.5rem) translateY(-0.5rem);
  transition: 0.5s;
}

.about__list-item:before {
  content: counter(list-counter);
  counter-increment: list-counter;
  width: 0.75rem;
  height: 0.75rem;
  padding: 0.5rem;
  margin-right: 0.5rem;
  border-radius: 1rem;
  border: 0.25em solid #ccc;
  background: #000;
  color: #fff;
  font-family: arial;
  font-weight: bold;
  text-align: center;
  display: inline-block;
}

.about__descbox {
  width: 70rem;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}

.about__descbox__column {
  width: 45%;
  height: 100%;
  margin: 0 1rem;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.about__descbox__column:last-child {
  align-items: center;
}

.about__desc-title {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 1.2rem;
  line-height: 1.4rem;
  color: rgb(50, 50, 50);
}

.about__desc-text {
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 0.9rem;
  line-height: 1.5rem;
  margin-top: 1rem;
  opacity: 0.7;
  color: rgb(40, 40, 40);
  margin-left: 0.4rem;
}

.about__desc-img {
  width: 100%;
  height: 12rem;
  object-fit: cover;
  border-radius: 1rem;
  margin: 1rem;
  margin-top: 3rem;
  margin-bottom: 3rem;
}

.about-section {
  width: 100%;
  padding: 2rem 0;
  padding-top: 5rem;
  display: flex;
  flex-direction: column;
  align-items: center;
  position: relative;
}

.about-section::before {
  content: "";
  background: linear-gradient(to bottom, #c3c3c3, #ffffff);
  width: 100%;
  height: 100%;
  position: absolute;
  left: 0;
  z-index: -1;
  top: 0;
}

.about-title-block {
  text-align: center;
  margin-bottom: 2rem;
}

.about-subtitle {
  font-size: 0.9rem;
  font-weight: 500;
  letter-spacing: 2px;
  color: #222;
  margin-bottom: 0.2rem;
}

.about-title {
  font-size: 2.5rem;
  font-weight: 800;
  color: #222;
  margin-bottom: 0.3rem;
}

.about-tagline {
  font-size: 1rem;
  font-style: italic;
  color: #555;
}

.about-container {
  width: 100%;
  display: flex;
  justify-content: center;
  margin-bottom: 2rem;
}

.about-content-box {
  background: #fff;
  padding: 2rem;
  border-radius: 1rem;
  max-width: 850px;
  box-shadow: 0 10px 20px rgba(0, 0, 0, 0.05);
  text-align: left;
}

.about-description p {
  font-size: 1rem;
  line-height: 1.8rem;
  color: #333;
  margin-bottom: 1.2rem;
}

.about-links-centered {
  display: flex;
  justify-content: center;
  gap: 1.5rem;
  flex-wrap: wrap;
  margin-top: 1.5rem;
}

.link-card {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-decoration: none;
  background: white;
  padding: 1rem;
  border-radius: 1rem;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08);
  transition: 0.3s ease;
  width: 140px;
}

.link-card:hover {
  transform: translateY(-5px);
  box-shadow: 0 8px 20px rgba(0, 0, 0, 0.15);
}

.link-card img {
  width: 48px;
  height: 48px;
  margin-bottom: 0.5rem;
  object-fit: contain;
}

.link-card span {
  font-weight: 600;
  color: #333;
  font-size: 0.95rem;
}

.sub-description__box {
  width: 100%;
  display: flex;
  justify-content: center;
  align-items: center;
}

.sub-description {
  width: 40%;
  color: #333333;
  /* 글자 색상 설정 */
}

.menu-list__largeCard {
  cursor: pointer;
  margin: 2rem 1rem;
  margin-bottom: 2rem;
  width: 26rem;
  height: 24rem;
  max-width: 80%;
  max-height: 80%;
  display: flex;
  align-items: center;
  justify-content: center;
  background: #f3f3f3d2;
  box-shadow: rgba(0, 0, 0, 0.34) 0px 0px 10px;
  padding: 1rem;
  border-radius: 1rem;
  /* transition: 0.5s; */
  position: relative;
}

.menu-list__largeCard:hover {
  /* width: 36rem;
  max-width: calc(80% + 2rem); */
  box-shadow: rgba(0, 0, 0, 0.31) 0px 5px 20px, rgba(0, 0, 0, 0.35) 0px 3px 10px;
  transition: 0.5s;
  background: #f3f3f3;
}

.menu-list__largeIcon {
  /* pointer-events: none; */
  width: 100%;
  height: 100%;
  object-fit: contain;
  position: absolute;
  opacity: 0.05;
  border-radius: 1rem;
  transition: opacity 1.5s ease;
  /* 투명도가 서서히 변화하는 트랜지션 설정 */
}

.menu-list__largeIcon:hover {
  opacity: 0.5;
}

.largeTextbox {
  width: 20rem;
  height: 6rem;
  text-align: center;
}

.menu-list__largeTitle {
  pointer-events: none;
  font-family: "Montserrat", sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 2rem;
  text-align: center;
  color: rgb(31, 56, 90);
}

/* Responsive Design */
@media (max-width: 1200px) {
  .menu-list__col {
    width: 45%;
  }
}

@media (max-width: 992px) {
  .menu-list__col {
    width: 100%;
  }

  .intro__descbox__column {
    flex: 1 1 100%;
  }

  .main__name {
    width: 80%;
    min-width: 20rem;
  }
}

@media (max-width: 768px) {
  .menu__title {
    font-size: 2rem;
  }

  .menu__description {
    font-size: 1.5rem;
    width: 90%;
  }

  .about-content-box {
    padding: 1.5rem;
  }

  .link-card {
    width: 120px;
  }
}

@media (max-width: 576px) {
  .main__name {
    width: 90%;
    min-width: 15rem;
  }

  .main__desc {
    font-size: 0.9rem;
  }

  .menu-list__card,
  .menu-list__largeCard {
    width: 90%;
  }

  .about-title {
    font-size: 2rem;
  }
}
</style>
