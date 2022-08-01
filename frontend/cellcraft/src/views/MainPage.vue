
<style scoped>
  @import '../../css/style.css';
  .info{
    display: flex;
    align-content: center;
    justify-content: center;
    align-items: center;
    background-color: #646FD4;
    color: white;
    padding-top: 4vh;
    flex-direction: column;
    padding-bottom: 4vh;
  }
  .info__name{
    font-size: 2.5vw;
  }
  .info__subname{
    margin-top: 3.5vh;
    font-size: .8vw;
    font-weight: bold;
  }
  .info__img img{
    margin-top: 3vh;
    background-color: white;
    width: 30vw;
    height: auto;
  }
  .infoTxt{
    background-color: #242F9B;
    color: white;
    text-align: center;
    padding: 2vh;
  }
  .menu{
    background-color: rgb(240,240,240);
    padding: 4vh 0;

  }
  .i-menu__txt{
    font-size: 1vw;
    text-align: center;
  }
  .i-menu__blocks{
    display: flex;
    justify-content: center;
    flex-wrap: wrap;
  }
  .menu__item{
    background: white;
    border: 0.1vw solid #C2C2C2;
    margin: 2vh 1.6vw;
    /* max-width: 15vw; */
    width: 15vw;
    height: auto;
    border-radius: .5vw;
    box-sizing: border-box;
    display: flex;
    padding: 2vh 1vw;
    align-items: center;
    word-break:keep-all;
  }
  .menu__item img{
    max-width: 3vw;
    max-height: 100%;
    margin-right: .5vw;
  }
  .menu__title{
    color: #269ABC;
    font-size: .9vw;
    font-weight: bold;
    margin-bottom: 1vh;
  }
  footer{
    background-color: #F8F9FA;
  }
  .copyright{
    padding: 2vh 0vh 1vw 1vw;
  }
  /* .tab{
    margin-top: 5vh;
  }
  .tab__menubar{
    border-bottom: solid .1vw #269ABC;
    display: flex;
    align-content: center;
    justify-content: center;
    margin: 0 10vw;
  }
  .tab__menu{
    color: #269ABC;
    padding: 0 1.8vw;
    padding-bottom: 1vh;
  }
  .tab__menu:hover, .tab__menu.scope{
    font-weight: bold;
    border-bottom: solid .8vw #269ABC;
  } */

</style>

<template>
  <!-- <div>
      <h1>Main</h1>
  </div> -->
  <section class="mainPage">
    <section class="info">
      <div class="info__name">
        Single Cell Analysis
      </div>
      <div class="info__subname">
        Inter-cell Type Gene Regulatory Network Model based on Integrative Analysis
      </div>
      <div class="info__img">
        <img src="@/assets/tenet.png" alt="">
      </div>
    </section>
    <section>
      <div class="infoTxt">
        Developed in the Data Sciences Platform at the Broad Institute, the toolkit offers a wide variety of tools with a primary focus on variant discovery and
        <br>
        genotyping. Its powerful processing engine and high-performance computing features make it capable of taking on projects of any size. Learn more
      </div>
    </section>
    <section class="menu">
      <div class="i-menu__txt">
        Find answers to your questions. Stay up to date on the latest topics. Ask questions and help others.
      </div>
      <ul class="i-menu__blocks">
        <li class="menu__item">
          <img src="@/assets/work-started.png" alt="">
          <div>
            <div class="menu__title">
              Getting Started
            </div>
            <div class="menu__info">
              Best practices, tutorials, and other info to get you started
            </div>
          </div>
        </li>
        <li class="menu__item">
          <img src="@/assets/documents-search.png" alt="">
          <div>
            <div class="menu__title">
              Technical Documentation
            </div>
            <div class="menu__info">
              Algorithms, glossary, and other detailed resources
            </div>
          </div>
        </li>
        <li class="menu__item">
          <img src="@/assets/loudspeaker.png" alt="">
          <div>
            <div class="menu__title">
              Announcements
            </div>
            <div class="menu__info">
              Blog and events
            </div>
          </div>
        </li>
      </ul>
      <ul class="i-menu__blocks">
        <li class="menu__item">
          <img src="@/assets/repair-fix-repairing.png" alt="">
          <div>
            <div class="menu__title">
              Tool Index
            </div>
            <div class="menu__info">
              Purpose, usage and options for each tool
            </div>
          </div>
        </li>
        <li class="menu__item">
          <img src="@/assets/faqs.png" alt="">
          <div>
            <div class="menu__title">
              Forum
            </div>
            <div class="menu__info">
              Ask our team for help and report issues
            </div>
          </div>
        </li>
        <li class="menu__item">
          <!-- <img src="../assets/stars.png" alt=""> -->
          <img src="@/assets/stars.png" alt="">
          <div>
            <div class="menu__title">
              Showcase
            </div>
            <div class="menu__info">
              Check out these fully configured workspaces
            </div>
          </div>
        </li>
      </ul>
    </section>
    <!-- <section class="tab">
      <div class="tab__menubar">
        <div class="tab__menu scope">
          Scope
        </div>
        <div class="tab__menu">
          Pipeline
        </div>
        <div class="tab__menu">
          Platforms
        </div>
        <div class="tab__menu">
          Tools/Code
        </div>
        <div class="tab__menu">
          License
        </div>
        <div class="tab__menu">
          Help
        </div>
      </div>
    </section> -->
    <footer>
      <div class="copyright">
       Copyright 2022. Cellcraft Co. all rights reserved
      </div>
    </footer>
  </section>
</template>

<script>
export default {
}
</script>

<style>

</style>
