<template>
  <div class="layout">
    <main class="main-component">
      <img class="main__img" src="@/assets/main_bg.png">
      <div class="main__textbox">
        <p class="main__name">CELLCRAFT</p>
        <span class="main__desc">
          Platform for Reconstruction of Inter-cell Type Gene Regulatory Network Model
          <br>
          based on Integrative Analysis of Single Cell RNA-seq, ATAC-seq and Spatial Transcriptomics Data
        </span>
      </div>
    </main>
    <div class="notice-component">
      <h1 class="notice__title">Notice</h1>
      <div class="notice__column">
        <p class="notice__info">[Notice] The lastest CELLCRAFT release came out...</p>
        <p class="notice__date">2023.01.04</p>
      </div>
      <div class="notice__column">
        <p class="notice__info">[Event] Welcome to CELLCRAFT   workshop</p>
        <p class="notice__date">2023.01.04</p>
      </div>
      <img class="notice__list" src="@/assets/menu.png">
    </div>
    <section class="menu-component">
      <ul class="menu-list">
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with CELLCRAFT</span>
          </div>
          <img class="arrow hidden" src="@/assets/right-arrow.png">
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/document.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Document</p>
            <span class="menu-list__desc">Algorithms, Data sets, Research results, detailed resources etc.</span>
          </div>
          <img class="arrow hidden" src="@/assets/right-arrow.png">
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/chat.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Forum</p>
            <span class="menu-list__desc">A community where users communicate freely with each other.</span>
          </div>
          <img class="arrow hidden" src="@/assets/right-arrow.png">
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/database.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Dataset</p>
            <span class="menu-list__desc">Single-cell data set accumulated through the study</span>
          </div>
          <img class="arrow hidden" src="@/assets/right-arrow.png">
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/megaphone.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Notice</p>
            <span class="menu-list__desc">Check out announcements, blogs, and events</span>
          </div>
          <img class="arrow hidden" src="@/assets/right-arrow.png">
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/faq.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">FAQ</p>
            <span class="menu-list__desc">You can see frequently asked questions while using our service</span>
          </div>
          <img class="arrow hidden" src="@/assets/right-arrow.png">
        </li>
      </ul>
    </section>
    <div class="intro-component">
      <div class="intro__textbox">
        <h1 class="intro__title">
          <span class="bold">ABOUT</span>
          <br>
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
          <h1 class="intro__desc-title">Single-cell analysis data platform for researchers</h1>
          <p class="intro__desc-text">
            CELLCRAFT is a data platform for easy analysis of single cell data and efficient use of data. CELLCRAFT contains many utilities for performing related tasks such as processing and quality control of high-throughput sequencing data and provides services on the web.
          </p>
          <p class="intro__desc-text">
            These tools are primarily designed for data analysts in bioinformatics, but even those who do not know programming can easily use CELLCRAFT to handle data analysis. CELLCRAFT was originally developed for single-cell analysis, but it is being developed to handle multiple biological data.
          </p>
        </div>
        <div class="intro__descbox__column">
          <img class="intro__desc-img" src="@/assets/example.jpg">
        </div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  methods: {
    mouseoverEvent (ev) {
      console.log(ev)
    },
    mouseleaveEvent (ev) {
      console.log(ev)
    }
  }
}
</script>

<style scoped>
.bold{
  font-weight: bold;
}
.arrow{
  display: inline-block;
  width: 1.5rem;
  height: 1.5rem;
  object-fit: cover;
  position: absolute;
  top: 1rem;
  right: 1rem;
}
.hidden{
  display: none;
}
.layout{
  width: 100%;
  height: 100%;
  position: relative;
  overflow: hidden;
}
.main-component,
.notice-component,
.menu-component,
.intro-component{
  width: 90rem;
  max-width: calc(100% - 3rem);
  margin-right: auto;
  margin-left: auto;
  position: relative;
}

.main-component{
  height: 20rem;
}
.notice-component{
  height: 5rem;
  margin-top: 2.5rem;
  display: flex;
  align-items: center;
}
.notice-component:before {
    content: '';
    background: #203969;
    width: 1000%;
    height: 100%;
    position: absolute;
    left: 0;
    margin-left: -500%;
    z-index: -1;
    top: 0;
}
.menu-component{
  padding: 2rem 0;
  display: flex;
  align-items: center;
  box-sizing: border-box;
}
.intro-component{
  height: 22.5rem;
  display: flex;
  align-items: center;
}
.intro-component:before {
  content: '';
  background: #f3f3f37a;
  width: 1000%;
  height: 100%;
  position: absolute;
  left: 0;
  margin-left: -500%;
  z-index: -1;
  top: 0;
}
.main__img{
  width: 90rem;
  max-width: calc(100%);
  height:inherit;
  object-fit: cover;
  object-position: left;
  border-radius: 1.5rem;
}
.main__textbox{
  width: 50rem;
  height: 8rem;
  position: absolute;
  top: calc(50% - 4rem);
  left: 5rem;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
}
.main__name{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 4rem;
  line-height: 4rem;
  text-align: left;
}
.main__desc{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1.04rem;
  line-height: 1.6rem;
  text-align: left;
}
.notice__title{
  width: 10%;
  height: 100%;
  display: flex;
  justify-content: flex-start;
  align-items: center;

  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 1.5rem;
  line-height: 1.5rem;
  text-align: left;
  color: #fff;
}
.notice__column{
  width: 45%;
  height: 100%;
  display: flex;
  align-items: center;
}
.notice__info,
.notice__date{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 1rem;
  line-height: 1rem;
  text-align: center;
  color: #fff;
  opacity: 0.7;
  display: flex;
  align-items: center;
  justify-content: center;
}
.notice__info{
  width: 70%;
  height: 100%;
}
.notice__date{
  width: 30%;
  height: 100%;
}
.notice__list{
  width: 1.5rem;
  height: 1.5rem;
  object-fit: cover;
  filter: invert(91%) sepia(91%) saturate(30%) hue-rotate(198deg) brightness(106%) contrast(100%);
}
.menu-list{
  width: 100%;
  padding: 1rem 0;
  overflow: hidden;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  justify-content: center;
}
.menu-list__card{
  cursor: pointer;
  margin: 0 1.5rem;
  margin-bottom: 2rem;
  width: 22.5rem;
  max-width: calc(33% - 5rem);
  height: 6.5rem;
  display: flex;
  align-items: center;
  /* box-shadow: rgba(0, 0, 0, 0.24) 0px 3px 8px; */
  box-shadow: rgba(0, 0, 0, 0.04) 0px 3px 5px;
  padding: 1rem;
  border-radius: 1rem;
  position: relative;
}
.menu-list__card:hover{
  box-shadow: rgba(0, 0, 0, 0.19) 0px 10px 20px, rgba(0, 0, 0, 0.23) 0px 6px 6px;
  transition: 0.5s;
}
.menu-list__icon{
  width: 3rem;
  height: 3rem;
  object-fit: cover;
  margin-right: 1rem;
}
.menu-list__textbox{
  width: 17.5rem;
  height: 5.5rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
  position: relative;
}
.menu-list__title{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 600;
  font-size: 1.4rem;
  line-height: 1.1rem;
  position: absolute;
  top: calc(100% - 5rem);
}
.menu-list__desc{
  width: calc(100%);
  max-width: 19rem;
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 0.9rem;
  line-height: 1.1rem;
  position: absolute;
  top: calc(100% - 3rem);
}
.intro__textbox{
  width: 20rem;
  height: 80%;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;
  border-right: 1px solid rgba(0, 0, 0, 0.3);
}
.intro__title{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 2.75rem;
  line-height: 2.75rem;
}
.intro__list{
  width: 100%;
  margin-top: 1rem;
  counter-reset: list-counter;
  list-style: none;
}
.intro__list-item{
  width: 70%;
  margin: 1rem 0;
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1rem;
  padding: 0.25rem;
}
.intro__list-item:hover{
  box-shadow: rgba(0, 0, 0, 0.4) -5px 5px, rgba(0, 0, 0, 0.3) -10px 10px;
  transform: translateX(0.5rem) translateY(-0.5rem);
  transition: 0.5s;
}
.intro__list-item:before{
  content: counter(list-counter);
  counter-increment: list-counter;
  width: 0.75rem;
  height: 0.75rem;
  padding: 0.5rem;
  margin-right: 0.5rem;
  border-radius: 50%;
  border: .25em solid #ccc;
  background: #000;
  color: #fff;
  font-family: arial;
  font-weight: bold;
  text-align: center;
  display: inline-block;
}
.intro__descbox{
  width: 70rem;
  height: 80%;
  display: flex;
  justify-content: center;
}
.intro__descbox__column{
  width: 45%;
  height: 100%;
  margin: 0 1rem;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
}
.intro__desc-title{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: bold;
  font-size: 1.2rem;
  line-height: 1.4rem;
}
.intro__desc-text{
  font-family: 'Montserrat', sans-serif;
  font-style: normal;
  font-weight: 400;
  font-size: 0.9rem;
  line-height: 1.5rem;
  margin-top: 1rem;
  opacity: 0.7;
}
.intro__desc-img{
  width: 80%;
  height: 100%;
  object-fit: cover;
  border-radius: 2rem;
}
</style>
