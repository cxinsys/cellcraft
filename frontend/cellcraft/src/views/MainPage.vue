<template>
  <div class="layout">
    <main class="main-component">
      <img class="main__img" src="@/assets/main_bg.png">
      <div class="main__textbox">
        <p class="main__name">CellCraft</p>
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
        <p class="notice__info">[Notice] The lastest CellCraft release came out...</p>
        <p class="notice__date">2023.01.04</p>
      </div>
      <div class="notice__column">
        <p class="notice__info">[Event] Wellcome to CellCraft workshop</p>
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
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
        <li class="menu-list__card">
          <img class="menu-list__icon" src="@/assets/e-book.png">
          <div class="menu-list__textbox">
            <p class="menu-list__title">Tutorial</p>
            <span class="menu-list__desc">The best tutorial you can see to get started with cellcraft</span>
          </div>
        </li>
      </ul>
    </section>
  </div>
</template>

<script>
export default {
}
</script>

<style scoped>
.layout{
  width: 100%;
  height: 100%;
}
.main-component,
.notice-component,
.menu-component{
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
.menu-component{
  padding: 2rem 0;
  display: flex;
  align-items: center;
  box-sizing: border-box;
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
  font-family: 'NanumGothic';
  font-style: normal;
  font-weight: bold;
  font-size: 4rem;
  line-height: 4rem;
  text-align: left;
}
.main__desc{
  font-family: 'NanumGothic';
  font-style: normal;
  font-weight: 500;
  font-size: 1.1rem;
  line-height: 1.4rem;
  text-align: left;
}
.notice__title{
  width: 10%;
  height: 100%;
  display: flex;
  justify-content: flex-start;
  align-items: center;
  font-family: 'NanumGothic';
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
  font-family: 'NanumGothic';
  font-style: normal;
  font-weight: 500;
  font-size: 1.1rem;
  line-height: 1.1rem;
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
}
.menu-list__icon{
  width: 3rem;
  height: 3rem;
  object-fit: cover;
  margin-right: 1rem;
}
.menu-list__textbox{
  height: 5.5rem;
  display: flex;
  flex-direction: column;
  justify-content: center;
}
.menu-list__title{
  font-family: 'NanumGothic';
  font-style: normal;
  font-weight: bold;
  font-size: 1.5rem;
  line-height: 1.1rem;
  margin-bottom: 1rem;
}
.menu-list__desc{
  font-family: 'NanumGothic';
  font-style: normal;
  font-weight: 500;
  font-size: 1rem;
  line-height: 1.1rem;
}
</style>
