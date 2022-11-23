
<style scoped>
  @import '../../css/style.css';
  .info{
    display: flex;
    align-content: center;
    justify-content: center;
    align-items: center;
    /* background-color: #646FD4; */
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
    width: 100vw;
    height: auto;
    background-color: rgb(240,240,240);
    /* padding-top: 4vh; */
    padding-bottom: 50px;
  }
  .i-menu__txt{
    font-size: 1vw;
    text-align: center;
  }
  .i-menu__blocks{
    display: flex;
    justify-content: center;
    flex-wrap: wrap;
    padding-top: 20px;
  }
  .menu__item{
    background: white;
    margin: 2vh 1.6vw;
    min-width: 130px;
    max-width: 240px;
    width: 20vw;
    max-height: 500px;
    min-height: 200px;
    height: 50vh;
    border-radius: .5vw;
    box-sizing: border-box;
    display: flex;
    padding: 2vh 1vw;
    align-items: center;
    word-break:keep-all;
    flex-direction: column;
    justify-content: space-evenly;
  }
  .menu__item:hover{
    box-shadow: 0 4px 5px rgba(0, 0, 0, 0.6);
  }
  .menu__item img{
    max-width: 50px;
    max-height: 100%;
    margin-right: .5vw;
    margin-bottom: 2vh;
  }
  .menu__title{
    color: #269ABC;
    font-size: 1.2rem;
    font-weight: bold;
    margin-bottom: 1vh;
    text-align: center;
  }
  .menu__info{
    margin-top: 20px;
  }
  .bg_img{
    background-color: rgba(0,0,0,0.7);
    display: flex;
    align-content: center;
    justify-content: center;
    align-items: center;
  }
  .figure_2 img {
    display: block;
    will-change: transform;
    transform: scale(calc(1 / var(--scale)));
    height: 60vh;
    width: 100%;
    max-width: 100%;
    min-width: 70%;
    object-fit: cover;
    object-position: center center;
    opacity: 0.55;
  }
  .wrapper{
    width:100%;
    position: relative;
  }
  .bg_txt{
    position: absolute;
    top: 45%;
    left: 50%;
    width: 100%;
    transform: translate(-50%, -50%);
    font-size: 3rem;
    text-align: center;
    z-index: 100;
    color: rgba(255, 255, 255, 1);
    /* color: rgba(1,1,1,1); */
    font-family: 'LABDigital';
    text-shadow : 3px 3px 6px black;
  }
  @font-face {
      font-family: 'LABDigital';
      font-weight: normal;
      font-style: normal;
      src: url('https://cdn.jsdelivr.net/gh/webfontworld/fontlab/LABDigital.eot');
      src: url('https://cdn.jsdelivr.net/gh/webfontworld/fontlab/LABDigital.eot?#iefix') format('embedded-opentype'),
          url('https://cdn.jsdelivr.net/gh/webfontworld/fontlab/LABDigital.woff2') format('woff2'),
          url('https://cdn.jsdelivr.net/gh/webfontworld/fontlab/LABDigital.woff') format('woff'),
          url('https://cdn.jsdelivr.net/gh/webfontworld/fontlab/LABDigital.ttf') format("truetype");
      font-display: swap;
  }
</style>

<template>
  <!-- <div>
      <h1>Main</h1>
  </div> -->
  <section class="mainPage" id="main">
      <figure class="figure_2">
        <div class="wrapper">
        <div class="bg_img">
          <img src="@/assets/main3.jpg" alt="" ref="ref_p">
        </div>
        <div class="bg_txt"><p>cell craft</p></div>
        </div>
      </figure>
          <section class="menu">
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
              <!-- <li class="menu__item">
                <img src="@/assets/loudspeaker.png" alt="">
                <div>
                  <div class="menu__title">
                    Forum
                  </div>
                  <div class="menu__info">
                    Ask our team for help and report issues
                  </div>
                </div>
              </li> -->
              <li class="menu__item" >
                <img src="@/assets/faqs.png" alt="">
                <div>
                  <div class="menu__title">
                    Forum
                  </div>
                  <div class="menu__info" >
                    Ask our team for help and report issues
                  </div>
                </div>
              </li>
            </ul>
            <!-- <ul class="i-menu__blocks">
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

                <img src="../assets/stars.png" alt=""> -->
                <!-- <img src="@/assets/stars.png" alt="">
                <div>
                  <div class="menu__title">
                    Showcase
                  </div>
                  <div class="menu__info">
                    Check out these fully configured workspaces
                  </div>
                </div>
              </li>
            </ul>-->

          </section>
    </section>
</template>

<script>
export default {
  data () {
    return {
      prevScrollTop: 0,
      nowScrollTop: 0,
      lastScrollTop: 0,
      img_fadeOut: false,
      img_width: '',
      fixed_img_width: ''
    }
  },
  methods: {
    wheelDelta () {
      return this.prevScrollTop - this.nowScrollTop > 0 ? 'up' : 'down'
    },
    scrollListener () {
      var tmp = document.documentElement.scrollTop
      this.nowScrollTop = document.documentElement.scrollTop
      if (this.wheelDelta() === 'down') {
        // console.log('down')
        // let val = this.$refs.ref_p.width
        // this.img_width = `${val}px`
        // this.fixed_img_width = `${this.img_width} - 3px`
        console.log(tmp)
        if (tmp > 1100 && tmp < 1250) {
          this.img_fadeOut = true
        }
      }
      if (this.wheelDelta() === 'up') {
        console.log(tmp)
        // if (tmp < 1200) {
        //   this.img_fadeOut = false
        // }
        // let val = this.$refs.ref_p.width
        // this.img_width = `${val}px`
        // this.fixed_img_width = `${this.img_width} + 3px`
        // this.img_fadeOut = false
      }

      this.prevScrollTop = this.nowScrollTop
      this.lastScrollTop = tmp
    }
  },
  created () {
    // window.addEventListener('scroll', this.scrollListener)
  }
}
</script>

<style>

</style>
