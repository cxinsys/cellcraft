<template>
  <div>
    <div class="nodeBox">
      <img class="nodeBox__icon" src="@/assets/heatMap2.png" />
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      title: "Result",
    };
  },
  methods: {},
};
</script>

<style scoped>
.nodeBox {
  width: 3rem;
  height: 3rem;
  display: flex;
  align-items: center;
  justify-content: center;
  text-align: center;
}
.nodeBox__icon {
  width: 70%;
  height: 70%;
  object-fit: contain;
  /* filter: invert(97%) sepia(99%) saturate(0%) hue-rotate(123deg) brightness(27%)
    contrast(101%); */
}

@media (prefers-color-scheme: dark) {
  /* .nodeBox__icon {
    filter: invert(97%) sepia(99%) saturate(0%) hue-rotate(123deg)
      brightness(107%) contrast(101%);
  } */
}
</style>
