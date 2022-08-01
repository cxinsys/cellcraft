<style scoped>
  .files__header{
    display: flex;
    align-items: center;
    margin-top: 50px;
    padding-bottom: 10px;
    border-bottom: 1px #a1a1a1 solid;
  }
  .files__header__name{
    display: flex;
    align-items: center;
    margin-left: 150px;
    font-size: 1vw;
    font-weight: bold;
  }
  .files__header__img{
    margin-left: 5px;
    width: 15px;
    height: 100%;

  }
  .files__header__size{
    margin-left: 600px;
    font-size: 1vw;
    font-weight: bold;
  }
  .files__header__date{
    margin-left: 300px;
    font-size: 1vw;
    font-weight: bold;
  }
</style>

<template>
  <section>
    <div class="files__header">
      <div class="files__header__name">
        <div> Name </div>
        <img class="files__header__img" src="@/assets/descending-filter-desc.png" alt="">
      </div>
      <div class="files__header__size"> Size </div>
      <div class="files__header__date"> Date </div>
    </div>
  </section>
</template>

<script>
export default {

}
</script>
